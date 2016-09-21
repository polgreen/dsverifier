/**
 * DSVerifier - Digital Systems Verifier (Quantization Error in Closed-loop)
 *
 * Federal University of Amazonas - UFAM
 *
 * Authors:       Hussama Ismail <hussamaismail@gmail.com>
 *
 * ------------------------------------------------------
 *
 * Verify the quantization error for digital systems in closed-loop.
 *
 * This property analyses the plant and the controller performance
 * when connected using SERIES or FEEDBACK. The verification engine
 * checks whether the digital controllers' FWL effects causes an unexpected
 * error percentual
 *
 * The engine consider nondet inputs and nondet zeroes states
 * for the desired realization (DFI, DFII, and TDFII).
 *
 * ------------------------------------------------------
*/

extern digital_system plant;
extern digital_system plant_cbmc;
extern digital_system controller;

int verify_error_closedloop(void){

	overflow_mode = WRAPAROUND;

	/* generating closed loop for series or feedback */
	float * c_num = controller.b;
	int c_num_size = controller.b_size;
	float * c_den = controller.a;
	int c_den_size = controller.a_size;

	/* quantizing controller coefficients */
	fxp_t c_num_fxp[controller.b_size];
	fxp_float_to_fxp_array(c_num, c_num_fxp, controller.b_size);
	fxp_t c_den_fxp[controller.a_size];
	fxp_float_to_fxp_array(c_den, c_den_fxp, controller.a_size);

	/* getting quantized controller coefficients  */
	float c_num_qtz[controller.b_size];
	fxp_to_float_array(c_num_qtz, c_num_fxp, controller.b_size);
	float c_den_qtz[controller.a_size];
	fxp_to_float_array(c_den_qtz, c_den_fxp, controller.a_size);

	/* getting plant coefficients */
	#if (BMC == ESBMC)
		float * p_num = plant.b;
		int p_num_size = plant.b_size;
		float * p_den = plant.a;
		int p_den_size = plant.a_size;
	#elif (BMC == CBMC)
		float * p_num = plant_cbmc.b;
		int p_num_size = plant.b_size;
		float * p_den = plant_cbmc.a;
		int p_den_size = plant.a_size;
	#endif

	float ans_num_float[100];
	float ans_num_qtz[100];
	int ans_num_size = controller.b_size + plant.b_size - 1;
	float ans_den_qtz[100];
	float ans_den_float[100];
	int ans_den_size = controller.a_size + plant.a_size - 1;

	#if (CONNECTION_MODE == SERIES)
		ft_closedloop_series(c_num_qtz, c_num_size, c_den_qtz, c_den_size, p_num, p_num_size, p_den, p_den_size, ans_num_qtz, ans_num_size, ans_den_qtz, ans_den_size);
		ft_closedloop_series(c_num, c_num_size, c_den, c_den_size, p_num, p_num_size, p_den, p_den_size, ans_num_float, ans_num_size, ans_den_float, ans_den_size);
	#elif (CONNECTION_MODE == FEEDBACK)
		ft_closedloop_feedback(c_num_qtz, c_num_size, c_den_qtz, c_den_size, p_num, p_num_size, p_den, p_den_size, ans_num_qtz, ans_num_size, ans_den_qtz, ans_den_size);
		ft_closedloop_feedback(c_num, c_num_size, c_den, c_den_size, p_num, p_num_size, p_den, p_den_size, ans_num_float, ans_num_size, ans_den_float, ans_den_size);
	#endif

	int i;
	float y_qtz[X_SIZE_VALUE];
	float y_float[X_SIZE_VALUE];
	float x_qtz[X_SIZE_VALUE];
	float x_float[X_SIZE_VALUE];
	float xaux_qtz[ans_num_size];
	float xaux_float[ans_num_size];

	/* prepare inputs (all possibles values in dynamical range) */
	float xaux[ans_num_size];
	float nondet_constant_input = nondet_float();
	__DSVERIFIER_assume(nondet_constant_input >= impl.min && nondet_constant_input <= impl.max);
	for (i = 0; i < X_SIZE_VALUE; ++i) {
		x_qtz[i] = nondet_constant_input;
		x_float[i] = nondet_constant_input;
		y_qtz[i] = 0;
		y_float[i] = 0;
	}
	for (i = 0; i < ans_num_size; ++i) {
		xaux_qtz[i] = nondet_constant_input;
		xaux_float[i] = nondet_constant_input;
	}

	float yaux_qtz[ans_den_size];
	float yaux_float[ans_den_size];
	float y0_qtz[ans_den_size];
	float y0_float[ans_den_size];

	int Nw = ans_den_size > ans_num_size ? ans_den_size : ans_num_size;
	float waux_qtz[Nw];
	float waux_float[Nw];
	float w0_qtz[Nw];
	float w0_float[Nw];

	#if (REALIZATION == DFI)
		for (i = 0; i < ans_den_size; ++i) {
			yaux_qtz[i] = 0;
			yaux_float[i] = 0;
		}
	#else
		for (i = 0; i < Nw; ++i) {
			waux_qtz[i] = 0;
			waux_float[i] = 0;
		}
	#endif

	float xk, temp;
	float *aptr, *bptr, *xptr, *yptr, *wptr;

	int j;
	for(i=0; i<X_SIZE_VALUE; ++i){

		/* direct form I realization */
		#if (REALIZATION == DFI)
			/* realization with controller quantized */
			shiftLDouble(x_qtz[i], xaux_qtz, ans_num_size);
			y_qtz[i] = float_direct_form_1(yaux_qtz, xaux_qtz, ans_den_qtz, ans_num_qtz, ans_den_size, ans_num_size);
			shiftLDouble(y_qtz[i], yaux_qtz, ans_den_size);
			/* realization with controller non quantized */
			shiftLDouble(x_float[i], xaux_float, ans_num_size);
			y_float[i] = float_direct_form_1(yaux_float, xaux_float, ans_den_float, ans_num_float, ans_den_size, ans_num_size);
			shiftLDouble(y_float[i], yaux_float, ans_den_size);
		#endif

		/* direct form II realization */
		#if (REALIZATION == DFII)
			/* realization with controller quantized */
			shiftRDfloat(0, waux_qtz, Nw);
			y_qtz[i] = float_direct_form_2(waux_qtz, x_qtz[i], ans_den_qtz, ans_num_qtz, ans_den_size, ans_num_size);
			/* realization with controller non quantized */
			shiftRDfloat(0, waux_float, Nw);
			y_float[i] = float_direct_form_2(waux_float, x_float[i], ans_den_float, ans_num_float, ans_den_size, ans_num_size);
		#endif

		/* transposed direct form II realization */
		#if (REALIZATION == TDFII)
		  /* realization with controller quantized */
			y_qtz[i] = float_transposed_direct_form_2(waux_qtz, x_qtz[i], ans_den_qtz, ans_num_qtz, ans_den_size, ans_num_size);
			/* realization with controller non quantized */
			y_float[i] = float_transposed_direct_form_2(waux_float, x_float[i], ans_den_float, ans_num_float, ans_den_size, ans_num_size);
		#endif

		/* error verification using a % setted by user */
		float __quant_error = ((fxp_to_float(y_qtz[i]) - y_float[i])/y_float[i]) * 100;
		__DSVERIFIER_assert(__quant_error < impl.max_error && __quant_error > (-impl.max_error));

	}

	return 0;
}
