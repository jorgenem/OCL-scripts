#!/bin/bash

f2py -h --overwrite-signature unfolding_test2_f2py.pyf -m unfolding_test2_f2py unfolding-test2-with_response_matrix_generation-20161123.f

f2py -c unfolding_test2_f2py.pyf unfolding-test2-with_response_matrix_generation-20161123.f