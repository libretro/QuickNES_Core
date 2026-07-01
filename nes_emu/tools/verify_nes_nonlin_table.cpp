/* Validation for the nonlinearizer DAC bake: recompute the curve from the
 * exact float formula and compare against the baked header. Must be identical.
 *   g++ -O2 -I.. -o v verify_nes_nonlin_table.cpp && ./v
 */
#include <cstdio>
#include <cstdint>
#include <cstring>
#include "nes_nonlin_table.h"

enum { table_bits = 11, table_size = 1 << table_bits };
static double nonlinear_tnd_gain(void) { return 0.75; }

int main(void)
{
	int16_t ref[table_size];
	float const gain = 0x7fff * 1.3f;
	int const range = (int) (table_size * nonlinear_tnd_gain());
	for (int i = 0; i < table_size; i++)
	{
		int const offset = table_size - range;
		int j = i - offset;
		float n = 202.0f / (range - 1) * j;
		float d = 0;
		if (n)
			d = gain * 163.67f / (24329.0f / n + 100.0f);
		ref[j & (table_size - 1)] = (int16_t)(int) d;
	}
	if (memcmp(ref, baked_nonlin_table, sizeof ref) == 0)
	{
		std::printf("nonlin bake: IDENTICAL (%d entries)\n", table_size);
		return 0;
	}
	std::printf("nonlin bake: MISMATCH\n");
	return 1;
}
