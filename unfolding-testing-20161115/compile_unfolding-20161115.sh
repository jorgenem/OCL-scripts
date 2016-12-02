#!/bin/bash

f2py -h --overwrite-signature unfolding_test1_f2py.pyf -m unfolding_test1_f2py unfolding-test1-20161115.f

f2py -c unfolding_test1_f2py.pyf unfolding-test1-20161115.f