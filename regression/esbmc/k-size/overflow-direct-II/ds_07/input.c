#include <dsverifier.h>

	digital_system ds = {
		.b = { 1.2670, -0.8493 },
		.b_size = 2,
		.a = { 1.0000000, -0.2543000 },
		.a_size = 2,
		.sample_time = 0.2
	};


		implementation impl = {
			.int_bits = 3,
			.frac_bits = 5,
			.min = -3.0,
			.max = 3.0
		};
