#makefile
test: 
	gcc md2.c lhep.c bn_ext.c  poly_eval.c error_hdl.c -I /usr/local/include/flint  -lgmp -lflint -lrelic -lrelic_s -o test
.PHONY : clean
 clean:
	rm -f test
	
