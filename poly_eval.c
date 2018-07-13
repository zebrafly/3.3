#include "md2.h"

// pka
typedef struct
{
  fq_poly_t B;
  fq_poly_t R;
} pk_a;

//eka
typedef struct
{
  int d;
  fq_poly_t A;
  fq_poly_t Q;
} ek_a;

//vkx
typedef struct
{
  g1_t B;
  g2_t R;
  g2_t pi;
} vk_x;

// this function generate a random polynomial that has small coefficients
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound)
{

  fq_poly_init(pol,ctx);

  bn_t z1;
  bn_null(z1);

  fmpz_t z2;
  fmpz_init(z2);

  fq_t z3;
  fq_init(z3,ctx);

  for (int n=0;n<deg;n++)
  {
      bn_rand_mod(z1,bound);
      bn2fmpz(z2,z1);
      fq_set_fmpz(z3,z2,ctx);
      fq_poly_set_coeff(pol,n,z3,ctx);
  }


  // release the memory
  bn_free(z1);
  fmpz_clear(z2);
  fq_clear(z3,ctx);
  return 0;

}

int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus)
{
   bn_mul_basic(out,in1,in2);
   bn_mod_basic(out,out,modulus);
   return 0;
}

// this function reduce the coefficients of a polynomial
int msg_modp(fq_poly_t m, lhe_par *par)
{

   long len=fq_poly_length(m,par->ctx);

   fq_t zq;
   fq_init(zq,par->ctx);

   fmpz_t zf;
   fmpz_init(zf);

   for (int i=0;i<len;i++)
   {

     fq_poly_get_coeff(zq,m,i,par->ctx);
     fmpz_set_str(zf,fq_get_str_pretty(zq,par->ctx), 10);
     fmpz_mod(zf,zf,par->pf);
     fq_set_fmpz(zq,zf,par->ctx);
     fq_poly_set_coeff(m,i,zq,par->ctx);
   }

    fq_clear(zq,par->ctx);
    fmpz_clear(zf);

   return 0;

}


// the vc_keygen algorithm
int vc_keygen(pk_a *pka, ek_a *eka, fq_poly_t A, lhe_par *par)
{
	// b0<-Fp*
	bn_t z1;
	bn_zero(z1);

	fmpz_t z2;
	fmpz_init(z2);

	fq_t b0;
	fq_init(b0,par->ctx);

	while (bn_is_zero(z1)) bn_rand_mod(z1,par->p);
	bn2fmpz(z2,z1);
    fq_set_fmpz(b0,z2,par->ctx);

	// B(X) = X^2+b0
	fq_t one;
	fq_init(one,par->ctx);
	fq_t zero;
	fq_init(zero,par->ctx);
	fq_poly_t B;
	fq_poly_init(B,par->ctx);
	fq_zero(zero,par->ctx);
	fq_one(one,par->ctx);
	fq_poly_set_coeff(B,2,one,par->ctx);
	fq_poly_set_coeff(B,1,zero,par->ctx);
	fq_poly_set_coeff(B,0,b0,par->ctx);
  // fq_poly_print_pretty(B,"x",par->ctx);
	// A = BQ+R
	fq_poly_t Q;
	fq_init(Q,par->ctx);
	fq_poly_t R;
	fq_init(R,par->ctx);
	fq_poly_divrem(Q,R,A,B,par->ctx);


	// pka
	fq_poly_init(pka->B,par->ctx);
	fq_poly_init(pka->R,par->ctx);
	fq_poly_set(pka->B,B,par->ctx);
	fq_poly_set(pka->R,R,par->ctx);

	// eka
	fq_poly_init(eka->A,par->ctx);
	fq_poly_init(eka->Q,par->ctx);
	fq_poly_set(eka->A,A,par->ctx);
	fq_poly_set(eka->Q,Q,par->ctx);

	return 0;
}

// the vc_pgen algorithm
int vc_pgen(pk_a *pka, vk_x *vkx, bn_t *x, lhe_par *par)
{
	fmpz_t z2;
	fmpz_init(z2);

	fq_t z3;
	fq_init(z3,par->ctx);

	bn2fmpz(z2,x);
    fq_set_fmpz(z3,z2,par->ctx);

	fq_t z;
  fq_init(z,par->ctx);
	// vkx->B = x^2+b0


	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,pka->B,z3,par->ctx);
	// fq_poly_print_pretty(pka->B,"x",par->ctx);
	// printf("\n\n");
  bn_t temp;
  fq2bn(temp,z,par->ctx);
	g1_mul(vkx->B,par->g,temp);
	g1_print(vkx->B);
	printf("\n\n");

	// vkx->R = r1*x+r0
	fq_init(z,par->ctx);
	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,pka->R,z3,par->ctx);
  fq2bn(temp,z,par->ctx);
	g2_mul(vkx->R,par->h,temp);
	g1_print(vkx->R);
	printf("\n\n");

	return 0;
}

// the vc_comp algorithm
int vc_comp(gt_t *y, vk_x *vkx, bn_t x, ek_a *eka, lhe_par *par)
{
	fmpz_t z2;
	fmpz_init(z2);

	fq_t z3;
	fq_init(z3,par->ctx);

	bn2fmpz(z2,x);
    fq_set_fmpz(z3,z2,par->ctx);

	fq_t z;


	// y = A(x) mod p

	fq_init(z,par->ctx);

	fq_poly_t m;
  fq_poly_init(m,par->ctx);
	fq_poly_set(m,eka->A,par->ctx);
  fq_poly_print_pretty(m,"x",par->ctx);
  printf("\n\n");
	msg_modp(m,par);
    fq_poly_print_pretty(m,"x",par->ctx);
    printf("\n\n");

	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,m,z3,par->ctx);
  bn_t temp;
  fq2bn(temp,z,par->ctx);
	gt_mul(y,par->gt,temp);
	gt_print(y);

	// vkx->pi = Q(x)
	g2_new(vkx->pi);
	fq_init(z,par->ctx);
	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,eka->Q,z3,par->ctx);
  fq2bn(temp,z,par->ctx);
	g2_mul(vkx->pi,par->g,temp);
	// fq_poly_print_pretty(eka->Q,"x",par->ctx);
	// g2_print(par->h);
	g1_print(vkx->pi);
	printf("\n\n");

	return 0;
}

// the vc_vrfy algorithm
int vc_vrfy(gt_t y, vk_x *vkx, lhe_par *par)
{
	// y?=B(x)*vkx->pi+R(x)
 	gt_t right;
 	gt_new(right);
    gt_set_unity(right);

    gt_mul(right,vkx->B,vkx->pi);
    g1_add(right,right,vkx->R);

	int temp =gt_cmp(y,right);

	if (temp == 2) return 2;
	return 0;
}


int main(int argc, char **argv)
{

// initialize the project /////////////////////////////////////////////////
        if (md2_init())
        {
                printf("Testing FAILED\n");
                printf("Problem initializing the library\n");
               return 1;
        }


  lhe_par *par=malloc(sizeof(*par));
  lhep_new(par);
printf("\n--------------------------------------begin test-------------------------\n\n\n");

        // rand polynomial algorithm
		fq_poly_t A;
		fq_poly_init(A,par->ctx);
        fq_poly_rand(A,64,par->ctx,par->p);


        // print



        // test the vc_keygen algorithm
        pk_a *pka=malloc(sizeof(*pka));
        ek_a *eka=malloc(sizeof(*eka));
        	fq_poly_init(pka->B,par->ctx);
          fq_poly_init(pka->R,par->ctx);
          fq_poly_init(eka->A,par->ctx);
          fq_poly_init(eka->Q,par->ctx);
        vc_keygen(pka,eka,A,par);
        // fq_poly_print_pretty(eka->A,"x",par->ctx);
        // printf("\n\n");
        // fq_poly_print_pretty(pka->B,"x",par->ctx);
        // printf("\n\n");
        // fq_poly_print_pretty(eka->Q,"x",par->ctx);
        // printf("\n\n");
        // fq_poly_print_pretty(pka->R,"x",par->ctx);

        // print



        // test the vc_pgen algorithm
        vk_x *vkx=malloc(sizeof(*vkx));
        	g1_new(vkx->B);
          	g2_new(vkx->R);
          	g2_new(vkx->pi);
        // generate a random message vector of d elements
        bn_t x;
        bn_rand_mod(x,par->p);

        // compute the problem instance
        vc_pgen(pka,vkx,x,par);

        // print


        // test the vc_comp algorithm
        gt_t y;
        gt_new(y);

        vc_comp(y,vkx,x,eka,par);
        // print

		// gt_print(y);



       // test the vc_vrfy algorithm
       int flag=vc_vrfy(y,vkx,par);
       printf("\n\n the verification result is flag=%d\n\n",flag);


printf("\n\n\n--------------------------------------end  test------------------------\n\n\n");


/*
        //release the memory
  bn_free(order);
  fq_clear(u,par->ctx);
  g2_free(zg);
  bn_free(zt);
  fq_poly_clear(s,par->ctx);
  fq_poly_clear(m,par->ctx);
  fq_poly_clear(m1,par->ctx);
  fq_poly_clear(m2,par->ctx);
*/


        return 0;
}
