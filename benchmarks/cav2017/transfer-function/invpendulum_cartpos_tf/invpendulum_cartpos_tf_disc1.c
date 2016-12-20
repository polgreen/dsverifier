#include <dsverifier.h>

implementation impl = { 
	.int_bits =  7,
	.frac_bits =   3,
	.max =  1.000000,
	.min =  -1.000000
	};

digital_system plant = { 
	.b = {  0 , 0.026311 , 0.01434 , -0.050254 , -0.038282 },
	.b_size =  5,
	.a = {  1 , -4.2553 , 6.5105 , -4.2553 , 1 },
	.a_size =  5 
	};

