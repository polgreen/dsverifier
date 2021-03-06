#include <dsverifier.h>

#define SCHEMA_COUNT 8

#if DS_ID == 1

	digital_system ds = {
		.b = { -2.70561932258292, 4.91891416097458, -2.98975424265897, 0.607457221322442 },
		.b_size = 4,
		.a = { 1.0, -0.246954439828182, -0.800014083515246, 0.356805858742342 },
		.a_size = 4,
		.sample_time = 0.5
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 2

	digital_system ds = {
		.b = { -10.6541798765226, 28.9820624061405, -26.2822459536125, 7.94546415262250 },
		.b_size = 4,
		.a = { 1.0, -0.976750653540241, -0.695199955503643, 0.688265939892838 },
		.a_size = 4,
		.sample_time = 0.1
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 3

	digital_system ds = {
		.b = { -19.2414736482723, 54.9701308152594, -52.3486145801513, 16.6178011176798 },
		.b_size = 4,
		.a = { 1.0, -1.20633950872401, -0.423011266643906, 0.633303983755908 },
		.a_size = 4,
		.sample_time = 0.05
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 4

	digital_system ds = {
		.b = { -28.8602335012427, 84.0786290054145, -81.6496235908929, 26.4305091092668 },
		.b_size = 4,
		.a = { 1.0, -1.39914676903959, -0.111714450559205, 0.512179344931653 },
		.a_size = 4,
		.sample_time = 0.03
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 5

	digital_system ds = {
		.b = { -58.8544737007610, 174.846413073101, -173.146281000398, 57.1542857242685 },
		.b_size = 4,
		.a = { 1.0, -1.92141324780054, 0.863510897744980, 0.0580048403364478 },
		.a_size = 4,
		.sample_time = 0.01
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 6

	digital_system ds = {
		.b = { -79.8561004850143, 238.400671160303, -237.238792985947, 78.6942127596070 },
		.b_size = 4,
		.a = { 1.0, -2.26995707001617, 1.54695565228205, -0.276981072005969 },
		.a_size = 4,
		.sample_time = 0.005
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 11,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 7

	digital_system ds = {
		.b = { -111.957457802305, 335.544346474481, -335.216643366135, 111.629754586207 },
		.b_size = 4,
		.a = { 1.0, -2.79568196431083, 2.59175994708314, -0.796077785225516 },
		.a_size = 4,
		.sample_time = 0.001
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 12,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 12,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 12,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif

#if DS_ID == 8

	digital_system ds = {
		.b = { -119.158714850888, 357.336454078930, -357.196818781631, 119.019079546242 },
		.b_size = 4,
		.a = { 1.0, -2.91304070183424, 2.82614887122274, -0.913108155920450 },
		.a_size = 4,
		.sample_time = 0.0004
	};

	#define	IMPLEMENTATION_COUNT 3

	#if IMPLEMENTATION_ID == 1
		implementation impl = {
			.int_bits = 12,
			.frac_bits = 4,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 2
		implementation impl = {
			.int_bits = 12,
			.frac_bits = 8,
			.min = -5.0,
			.max = 5.0
		};
	#endif

	#if IMPLEMENTATION_ID == 3
		implementation impl = {
			.int_bits = 12,
			.frac_bits = 12,
			.min = -5.0,
			.max = 5.0
		};
	#endif

#endif
