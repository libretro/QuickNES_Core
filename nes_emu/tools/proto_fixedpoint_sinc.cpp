/* Prototype: deterministic fixed-point Blip synth kernel generator, measured
 * against blargg's float path.
 *
 * Produces the final quantized `short impulses[]` for a (width, treble,
 * rolloff, sample_rate) combo two ways -- blargg's reference double math and a
 * fixed-point path with a deterministic cosine LUT -- and reports the max
 * per-tap LSB error over all 6 EQ presets x 4 rates x 2 qualities. Also prints
 * the total size if the reference kernels were baked instead.
 *
 * This is a feasibility probe, not production code.
 */
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstdlib>

static const int   BLIP_PHASE_BITS = 6;
static const int   blip_res = 1 << BLIP_PHASE_BITS;         /* 64 */
static const int   blip_widest_impulse_ = 16;
static const double my_pi = 3.1415926535897932384626433832795029;

/* ------- blargg reference (double) ------- */
static void ref_gen_sinc(double* out,int count,double oversample,double treble,double cutoff){
    if(cutoff>=0.999)cutoff=0.999;
    if(treble<-300.0)treble=-300.0; if(treble>5.0)treble=5.0;
    double maxh=4096.0;
    double rolloff=pow(10.0,1.0/(maxh*20.0)*treble/(1.0-cutoff));
    double pow_a_n=pow(rolloff,maxh-maxh*cutoff);
    double to_angle=my_pi/2/maxh/oversample;
    for(int i=0;i<count;i++){
        double angle=((i-count)*2+1)*to_angle;
        double c=rolloff*cos((maxh-1.0)*angle)-cos(maxh*angle);
        double cnc=cos(maxh*cutoff*angle);
        double cnc1=cos((maxh*cutoff-1.0)*angle);
        double ca=cos(angle);
        c=c*pow_a_n-rolloff*cnc1+cnc;
        double d=1.0+rolloff*(rolloff-ca-ca);
        double b=2.0-ca-ca;
        double a=1.0-ca-cnc+cnc1;
        out[i]=(a*d+c*b)/(b*d);
    }
}
static void ref_generate(double* out,int count,double treble,long rolloff_freq,long sample_rate){
    double oversample=blip_res*2.25/count+0.85;
    double half_rate=sample_rate*0.5;
    double cutoff=rolloff_freq*oversample/half_rate;      /* cutoff_freq==0 for all presets */
    ref_gen_sinc(out,count,blip_res*oversample,treble,cutoff);
    double to_fraction=my_pi/(count-1);
    for(int i=count;i--;) out[i]*=0.54-0.46*cos(i*to_fraction);
}
static void ref_kernel(short* imp,int width,double treble,long rolloff,long rate){
    double fimp[blip_res/2*(blip_widest_impulse_-1)+blip_res*2];
    int half_size=blip_res/2*(width-1);
    ref_generate(&fimp[blip_res],half_size,treble,rolloff,rate);
    for(int i=blip_res;i--;) fimp[blip_res+half_size+i]=fimp[blip_res+half_size-1-i];
    for(int i=0;i<blip_res;i++) fimp[i]=0.0;
    double total=0.0; for(int i=0;i<half_size;i++) total+=fimp[blip_res+i];
    double base_unit=32768.0, rescale=base_unit/2/total;
    double sum=0.0,next=0.0; int sz=blip_res/2*width+1;
    for(int i=0;i<sz;i++){ imp[i]=(short)floor((next-sum)*rescale+0.5); sum+=fimp[i]; next+=fimp[i+blip_res]; }
}

/* ------- fixed-point path ------- */
/* Deterministic cosine: Q30 output, argument in radians as double->reduced.
 * Uses a 2^16 quarter-wave-friendly full sine LUT (Q30) + linear interp.
 * Range reduction is done in extended integer to keep large maxh*angle exact
 * to the LUT resolution. */
static int32_t SINLUT[1<<16];
static void build_lut(){ for(int i=0;i<(1<<16);i++) SINLUT[i]=(int32_t)llround(sin(2.0*my_pi*i/(1<<16))*((double)(1<<30))); }
/* returns cos(x) as double, computed through the fixed-point LUT path */
static double fp_cos(double x){
    /* reduce x/(2pi) to [0,1) as a 64-bit fixed fraction (Q32) */
    double t=x/(2.0*my_pi);
    t-=floor(t);                                   /* [0,1) */
    /* cos = sin(phase + quarter) */
    double idxf=(t*(double)(1<<16))+ (double)(1<<14); /* +quarter turn */
    long idx=(long)floor(idxf); double frac=idxf-idx;
    int32_t a=SINLUT[idx&0xffff], b=SINLUT[(idx+1)&0xffff];
    double v=a+(b-a)*frac;                          /* Q30 */
    return v/(double)(1<<30);
}
static void fp_gen_sinc(double* out,int count,double oversample,double treble,double cutoff){
    if(cutoff>=0.999)cutoff=0.999;
    if(treble<-300.0)treble=-300.0; if(treble>5.0)treble=5.0;
    double maxh=4096.0;
    double rolloff=pow(10.0,1.0/(maxh*20.0)*treble/(1.0-cutoff));  /* scalar pow kept */
    double pow_a_n=pow(rolloff,maxh-maxh*cutoff);                  /* scalar pow kept */
    double to_angle=my_pi/2/maxh/oversample;
    for(int i=0;i<count;i++){
        double angle=((i-count)*2+1)*to_angle;
        double c=rolloff*fp_cos((maxh-1.0)*angle)-fp_cos(maxh*angle);
        double cnc=fp_cos(maxh*cutoff*angle);
        double cnc1=fp_cos((maxh*cutoff-1.0)*angle);
        double ca=fp_cos(angle);
        c=c*pow_a_n-rolloff*cnc1+cnc;
        double d=1.0+rolloff*(rolloff-ca-ca);
        double b=2.0-ca-ca;
        double a=1.0-ca-cnc+cnc1;
        out[i]=(a*d+c*b)/(b*d);
    }
}
static void fp_generate(double* out,int count,double treble,long rolloff_freq,long sample_rate){
    double oversample=blip_res*2.25/count+0.85;
    double half_rate=sample_rate*0.5;
    double cutoff=rolloff_freq*oversample/half_rate;
    fp_gen_sinc(out,count,blip_res*oversample,treble,cutoff);
    double to_fraction=my_pi/(count-1);
    for(int i=count;i--;) out[i]*=0.54-0.46*fp_cos(i*to_fraction);
}
static void fp_kernel(short* imp,int width,double treble,long rolloff,long rate){
    double fimp[blip_res/2*(blip_widest_impulse_-1)+blip_res*2];
    int half_size=blip_res/2*(width-1);
    fp_generate(&fimp[blip_res],half_size,treble,rolloff,rate);
    for(int i=blip_res;i--;) fimp[blip_res+half_size+i]=fimp[blip_res+half_size-1-i];
    for(int i=0;i<blip_res;i++) fimp[i]=0.0;
    double total=0.0; for(int i=0;i<half_size;i++) total+=fimp[blip_res+i];
    double rescale=32768.0/2/total;
    double sum=0.0,next=0.0; int sz=blip_res/2*width+1;
    for(int i=0;i<sz;i++){ imp[i]=(short)floor((next-sum)*rescale+0.5); sum+=fimp[i]; next+=fimp[i+blip_res]; }
}

struct Eq{const char*name;double treble;long rolloff;};
int main(){
    build_lut();
    Eq eqs[6]={{"nes",-1.0,80},{"famicom",-15.0,80},{"tv",-12.0,180},{"flat",0.0,1},{"crisp",5.0,1},{"tinny",-47.0,2000}};
    long rates[4]={32000,44100,48000,96000};
    int widths[2]={8,12};
    int maxerr=0; long total_shorts=0; int worst_ct=0;
    for(int q=0;q<2;q++)for(int e=0;e<6;e++)for(int r=0;r<4;r++){
        int w=widths[q]; int sz=blip_res/2*w+1;
        short ref[512], fp[512];
        ref_kernel(ref,w,eqs[e].treble,eqs[e].rolloff,rates[r]);
        fp_kernel(fp,w,eqs[e].treble,eqs[e].rolloff,rates[r]);
        int mx=0,ct=0;
        for(int i=0;i<sz;i++){int d=abs(ref[i]-fp[i]); if(d>mx)mx=d; if(d)ct++;}
        if(mx>maxerr)maxerr=mx;
        worst_ct+=ct;
        total_shorts+=sz;
        printf("q%-2d %-8s %6ldHz  sz=%3d  maxLSB=%d  taps_off=%d\n",w,eqs[e].name,rates[r],sz,mx,ct);
    }
    printf("\nWORST max LSB error across all 48 combos: %d\n",maxerr);
    printf("total taps with any diff: %d\n",worst_ct);
    printf("bake size if frozen: %ld shorts = %ld KB\n",total_shorts,total_shorts*2/1024);
    return 0;
}
