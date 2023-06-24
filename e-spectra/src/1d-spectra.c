
#include "definitions.h"
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	int i,j,k,a;
	int xd,yd,zd;
	int nproc,rank;
	PDATA pd;
	POINTER V;
	int tst;
	
	char fn[80];
	FILE *fd;
	MPI_File fh;
	MPI_Offset offset;
	MPI_Status status;
	
	ptrdiff_t N0, N1;
	fftw_plan plan,inv_plan,plan_u,plan_v,plan_w;
	fftw_complex *data;
	double *rin_u,*rin_v,*rin_w;
	fftw_complex *cout_u,*cout_v,*cout_w;
	ptrdiff_t alloc_local, local_n0, local_n0_start, local_n1_aft, local_n1_start_aft;
	
	double euux,euuy,evvx,evvy,ewwx,ewwy;
	double *Euux,*Euuy,*Evvx,*Evvy,*Ewwx,*Ewwy;
	double *tempy;
	double fac,fac2;
	
	MPI_Init(&argc, &argv);
	fftw_mpi_init();
	
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	pd.myrank=rank;
	pd.numproc=nproc;
	
	if (pd.myrank != 0)
		pd.left = pd.myrank - 1;
	else
		pd.left = pd.numproc - 1;
	
	if (pd.myrank != (pd.numproc-1))
		pd.right = pd.myrank + 1;
	else
		pd.right = 0;

	zd=zdv+2;
	yd=ydv;
	xd=1+(xdv/pd.numproc);
	
	N0 = xdv;
	N1 = ydv;
	
	local_n0 = N0/nproc;
	local_n0_start = rank*local_n0;

	local_n1_aft = (N1/2 +1)/nproc;
	local_n1_start_aft = rank*local_n1_aft;
	
	/* Allocating the velocity field */
	V.vel = (DP *)calloc((xd-1)*yd*zd*3,sizeof(DP));
	
	/* Allocating the E variables */
	Euux = (DP *)calloc(local_n0*zd,sizeof(DP));
	Evvx = (DP *)calloc(local_n0*zd,sizeof(DP));
	Ewwx = (DP *)calloc(local_n0*zd,sizeof(DP));
	
	Euuy = (DP *)calloc((N1/2 +1)*zd,sizeof(DP));
	Evvy = (DP *)calloc((N1/2 +1)*zd,sizeof(DP));
	Ewwy = (DP *)calloc((N1/2 +1)*zd,sizeof(DP));
	
	tempy = (DP *)calloc((N1/2 +1)*zd,sizeof(DP));
  
	/* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_2d(N0, (N1/2 +1), MPI_COMM_WORLD, &local_n0, &local_n0_start);
	
	rin_u = fftw_alloc_real(2 * alloc_local);
	rin_v = fftw_alloc_real(2 * alloc_local);
	rin_w = fftw_alloc_real(2 * alloc_local);
	
	cout_u = fftw_alloc_complex(alloc_local);
	cout_v = fftw_alloc_complex(alloc_local);
	cout_w = fftw_alloc_complex(alloc_local);
     
        /* create plan for in-place forward DFT */
	plan_u = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_u, cout_u, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_v = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_v, cout_v, MPI_COMM_WORLD, FFTW_ESTIMATE);
	plan_w = fftw_mpi_plan_dft_r2c_2d(N0, N1, rin_w, cout_w, MPI_COMM_WORLD, FFTW_ESTIMATE);
//	inv_plan = fftw_mpi_plan_dft_c2r_2d(N0, N1, cout,rin, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
	
	/* main loop */
	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V = reading(xd,yd,zd,V,pd,tst);
		
		for (k=0;k<zd;k++)
		{
			for (i=0;i<local_n0;i++)
			{
				Euux[i*zd+k] = 0.;
				Evvx[i*zd+k] = 0.;
				Ewwx[i*zd+k] = 0.;
			}
			for (j=0;j<(N1/2+1);j++)
			{
				Euuy[j*zd+k] = 0.;
				Evvy[j*zd+k] = 0.;
				Ewwy[j*zd+k] = 0.;
			}
		}
		
		/* moving the data into the fft plane */
		for (k=1;k<(zd-1);k++)
		{
			for (i=0;i<local_n0;i++)
			{
				for (j=0;j<N1;j++)
				{
					rin_u[i*(2*(N1/2 +1))+j] = Fb(V.vel,i,j,k,0,xd,yd,zd,3);
					rin_v[i*(2*(N1/2 +1))+j] = Fb(V.vel,i,j,k,1,xd,yd,zd,3);
					rin_w[i*(2*(N1/2 +1))+j] = Fb(V.vel,i,j,k,2,xd,yd,zd,3);
				}
			}
			/* compute transforms, in-place, as many times as desired */
			fftw_execute(plan_u);
			fftw_execute(plan_v);
			fftw_execute(plan_w);
			
			/* Scale the R => F FFT */
			for (i=0; i<local_n0; i++)
			{
				for (j=0;j<(N1/2 +1); j++)
				{
					cout_u[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_u[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_v[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_v[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
					
					cout_w[i*(N1/2 +1) +j][0] *= 2.*(1./((double)(2.*N0*N1)));
					cout_w[i*(N1/2 +1) +j][1] *= 2.*(1./((double)(2.*N0*N1)));
				}
			}
			
			/* subtract the average */
			if (pd.myrank == 0)
			{
				cout_u[0][0] = 0.; cout_u[0][1] = 0.;
				cout_v[0][0] = 0.; cout_v[0][1] = 0.;
				cout_w[0][0] = 0.; cout_w[0][1] = 0.;
			}
			
			/* Calculation of the energy spectra */
			/* Y spectra */
			for (j=0;j<(N1/2 + 1 );j++)
			{				
				for (i=0; i< local_n0 ;i++)
				{
					if (j==0)
						fac = 1.;
					else
						fac = 2.;

					Euux[i*zd+k] += fac*0.5*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvx[i*zd+k] += fac*0.5*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwx[i*zd+k] += fac*0.5*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
					
					Euuy[j*zd+k] += 0.5*fac*(cout_u[i*(N1/2 +1) +j][0]*cout_u[i*(N1/2 +1) +j][0] + cout_u[i*(N1/2 +1) +j][1]*cout_u[i*(N1/2 +1) +j][1]);
					Evvy[j*zd+k] += 0.5*fac*(cout_v[i*(N1/2 +1) +j][0]*cout_v[i*(N1/2 +1) +j][0] + cout_v[i*(N1/2 +1) +j][1]*cout_v[i*(N1/2 +1) +j][1]);
					Ewwy[j*zd+k] += 0.5*fac*(cout_w[i*(N1/2 +1) +j][0]*cout_w[i*(N1/2 +1) +j][0] + cout_w[i*(N1/2 +1) +j][1]*cout_w[i*(N1/2 +1) +j][1]);
				}
			  
			}
		}
			
		/* Mixing the y direction of all domains */
		MPI_Reduce(Euuy,tempy,(zd*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Euuy,tempy,(zd*(N1/2+1))*sizeof(double));
		MPI_Reduce(Evvy,tempy,(zd*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Evvy,tempy,(zd*(N1/2+1))*sizeof(double));
		MPI_Reduce(Ewwy,tempy,(zd*(N1/2+1)),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		memcpy(Ewwy,tempy,(zd*(N1/2+1))*sizeof(double));
		
		/* Exporting the y results: only by rank 0 */
		if (pd.myrank == 0)
		{
			sprintf(fn,"Euuy.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Euuy,sizeof(double),zd*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Evvy.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Evvy,sizeof(double),zd*(N1/2+1),fd);
			fclose(fd);
			
			sprintf(fn,"Ewwy.%.4d",tst);
			fd=fopen(fn,"wb");
			fwrite(Ewwy,sizeof(double),zd*(N1/2+1),fd);
			fclose(fd);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
			
		/* Exporting the x results, by all in parallel */
		sprintf(fn,"Euux.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zd),Euux,(local_n0*zd),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Evvx.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zd),Evvx,(local_n0*zd),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
		sprintf(fn,"Ewwx.%.4d",tst);
		MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
		MPI_File_set_view(fh,0,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
		MPI_File_write_at_all(fh,(pd.myrank*local_n0*zd),Ewwx,(local_n0*zd),MPI_DOUBLE,&status);
		MPI_File_close(&fh);
		
	}
     

     
        
	/* Destroying the plans */
	fftw_destroy_plan(plan_u);
	fftw_destroy_plan(plan_v);
	fftw_destroy_plan(plan_w);
//	fftw_destroy_plan(inv_plan);
	
	free(Euux);
	free(Euuy);
	free(Evvx);
	free(Evvy);
	free(Ewwx);
	free(Ewwy);
	free(tempy);
	
	MPI_Finalize();
  
 return 0; 
}
