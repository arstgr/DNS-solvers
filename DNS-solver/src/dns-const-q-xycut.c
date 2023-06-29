/***************************************************************************
 *   Amirreza Rastegari                                                    *
 *   arstgr@gmail.com                                                      *
 *                                                                         *
 *   Parallel LBM-BGK code with constant flow rate Q for DNS of turbulent  *
 *   channel flows                                                         *
 *   This program uses a new formulation for the constant flow rate        *
 *   implementation                                                        *
 ***************************************************************************/
#include "definitions.h"
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	
	int xd,yd,zd,TM;
	int i,j,k,q,a,b,c,nut;
	int tagl=1,tagr=2;
	int cntr,code;
	long ss,up;
	PDATA pd;
	FILE *epp;
	POINTER V;
	const DP C=1.;
	DP dx=1.;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	DP tau=0.;
	DP rhozero=1.,rhozeroinv;
	DP t1,t2,t3,t4,t5;
	char ch;
	FILE *mysave,*tm;
	FILE *sv,*pg,*vfx,*vfy,*vfz;
	char fn[60],fn2[40];
	DP utime,umax,dt;
	DP ux,uy,uz,dens;
	DP Ret, Reb, ut,dplus,Lplus,ubulk;

	/* MPI variables */
	MPI_File fh;
	MPI_Status status;
	MPI_Offset offset,offsett;
	MPI_Comm cart_grid;
	int dim[2], period[2], reorder;
	int coord[2], neigh[2], id;
	int rank, size;
	MPI_Datatype distf,block,sarray,larray;
	MPI_Datatype send_to_dw,recv_from_up, send_to_le,recv_from_ri,send_to_dl, recv_from_ur;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2],dimension;
	MPI_Info info;

	PI=4.*atan(1.);
	
	/* New Way of Computing input parameters (always add with 2) */
	zd=zdv+2;

	/* defining the inverse of rhozer: for performance issues */
	rhozeroinv=1./rhozero;
	
	/* MPI Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", STRIPE_COUNT);
	MPI_Info_set(info, "cb_nodes", CB_NODES_FCT);
//  MPI_Info_set(info, "striping_unit", STRIPE_SIZE);
    MPI_Info_set(info, "romio_cb_write", "enable");
//  MPI_Info_set(info, "romio_ds_write", "disable");
	
	pd.myrank = rank;
	pd.numproc = size;
	
	if (XDIM * YDIM != pd.numproc)
	{
		printf("Incorrect dimensions for the topology\n");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	dim[0]=XDIM; dim[1]=YDIM;
	period[0]=1; period[1]=1;
	reorder=1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &cart_grid);
	MPI_Cart_coords(cart_grid, pd.myrank, 2, coord);
	
	neigh[0] = coord[0];
	neigh[1] = coord[1] + 1;
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.up = rank;
	
	neigh[0] = coord[0];
	neigh[1] = coord[1] - 1;
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.dw = rank;
	
	neigh[0] = coord[0] + 1;
	neigh[1] = coord[1];
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.right = rank;
	
	neigh[0] = coord[0] - 1;
	neigh[1] = coord[1];
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.left = rank;
	
	neigh[0] = coord[0] - 1;
	neigh[1] = coord[1] - 1;
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.dl = rank;
	
	neigh[0] = coord[0] + 1;
	neigh[1] = coord[1] + 1;
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.ur = rank;
	//
	neigh[0] = coord[0] - 1;
	neigh[1] = coord[1] + 1;
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.ul = rank;
	
	neigh[0] = coord[0] + 1;
	neigh[1] = coord[1] - 1;
	MPI_Cart_rank(cart_grid, coord, &rank);
	pd.dr = rank;
	
	MPI_Type_contiguous(19,MPI_DOUBLE,&distf);
	MPI_Type_contiguous(zd,distf,&block);
	
	MPI_Type_commit(&block);
	
	/* no. of rows and columns in global array*/ 
	gsizes[0] = xdv;    gsizes[1] = ydv;     
	/* no. of processors in x and y directions */ 
	psizes[0] = XDIM;    psizes[1] = YDIM;
	/* size of the local array without the ghost cells */
	lsizes[0] = xdv/XDIM; lsizes[1] = ydv/YDIM;
	/* start of the array in terms of global array indices, ignoring the ghost cells */
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	
	/* Creating the type for the arrays withou ghost cells
	 * the idea is to read them with a type of array with ghost cells */
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, block, &sarray);
	MPI_Type_commit(&sarray);
	
	/* size of the local array with the ghost cells */
	memsizes[0] = lsizes[0] + 1;
	memsizes[1] = lsizes[1] + 1;
	
	/* start of the local array with ghost cells: note: ghost cells are the last row and columns */
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &larray);
	MPI_Type_commit(&larray);
	
	/* data type for transfer in y direction */
	lsizes[0] = xdv/XDIM ;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dw);
	MPI_Type_commit(&send_to_dw);
	
	array_start[0] = 0;
	array_start[1] = (ydv/YDIM);
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_from_up);
	MPI_Type_commit(&recv_from_up);
	
	/* data type for transfer in x direction */
	lsizes[0] = 1;
	lsizes[1] = ydv/YDIM;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_le);
	MPI_Type_commit(&send_to_le);
	
	array_start[0] = xdv/XDIM;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_from_ri);
	MPI_Type_commit(&recv_from_ri);
	
	/* data type for transfer to the corners */
	lsizes[0] = 1;
	lsizes[1] = 1;
	array_start[0] = 0;
	array_start[1] = 0;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &send_to_dl);
	MPI_Type_commit(&send_to_dl);
	
	array_start[0] = xdv/XDIM;
	array_start[1] = ydv/YDIM;
	MPI_Type_create_subarray(2, memsizes, lsizes, array_start, MPI_ORDER_C, block, &recv_from_ur);
	MPI_Type_commit(&recv_from_ur);

	t1=MPI_Wtime();

	/* Computing solution parameters */

	Reb=Rebs;
	ubulk=ub;

	yd=1+(ydv/YDIM);
	
	xd=1+(xdv/XDIM);

	/* Computing the Maximum velocity in lattice */
	umax=1.;//ut*(2.5*log(Ret)+5.5);

	/* Computing the required iterations for one non-dimensional time */
	utime=((double)(0.5*(zd-2.)/CFL));
	dt=1./utime;
	utime=(int)utime;

	i=xd/2;j=yd/2;k=zd/2;
	tau=(3./Reb)*ubulk*(zd-2.)+0.5;
	tau=1./tau;
	Fx=0.;
	Fy=0.;
	Fz=0.;

	/* defining the effective parameters in collision for performance issues */
	Gx=Fx;
	Gy=0.;
	Gz=0.;
	
	offsett = (MPI_Offset)((xd-1)*(yd-1)*zd);
	offsett *= ((MPI_Offset)19);
	offsett *= ((MPI_Offset)sizeof(double));

	/*Maximum loop iterations */
		nut=TEND-TSTR;
		up=nut*utime;

#ifdef ITER
		up=ITER;
#endif

	/* Memory Alocation */
	V.s=(int *)calloc(xd*yd*zd,sizeof(int));
	if (V.s == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.s: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.f=(DP *)calloc(xd*yd*zd*19,sizeof(DP));
	if (V.f == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.f: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	V.ftemp=(DP *)calloc(xd*yd*zd*19,sizeof(DP));
	if (V.ftemp == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.ftemp: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.rightbufs=(DP *)calloc(5*yd*zd,sizeof(DP));
	if (V.rightbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.rightbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.rightbufr=(DP *)calloc(5*yd*zd,sizeof(DP));
	if (V.rightbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.rightbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.leftbufs=(DP *)calloc(5*yd*zd,sizeof(DP));
	if (V.leftbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.leftbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.leftbufr=(DP *)calloc(5*yd*zd,sizeof(DP));
	if (V.leftbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.leftbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.upbufs=(DP *)calloc(5*xd*zd,sizeof(DP));
	if (V.upbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.upbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.upbufr=(DP *)calloc(5*xd*zd,sizeof(DP));
	if (V.upbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.upbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.dwbufs=(DP *)calloc(5*xd*zd,sizeof(DP));
	if (V.dwbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dwbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.dwbufr=(DP *)calloc(5*xd*zd,sizeof(DP));
	if (V.dwbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dwbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	
	V.urbufs=(DP *)calloc(9*zd,sizeof(DP));
	if (V.urbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.urbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.dlbufr=(DP *)calloc(9*zd,sizeof(DP));
	if (V.dlbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dlbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.dlbufs=(DP *)calloc(9*zd,sizeof(DP));
	if (V.dlbufs == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.dlbufs: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}
	V.urbufr=(DP *)calloc(9*zd,sizeof(DP));
	if (V.urbufr == NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.urbufr: %s\n",strerror(errno));
		MPI_Finalize();
		return -1;
	}

	for (i=0;i<12;i++)
		V.req[i] = MPI_REQUEST_NULL;
	/*Two first req are for send (left-right) and the second two are for recv (left-right)
	so even numbers are for right and odd for left */

	/* Establishing a channel for communication */
	MPI_Send_init(V.leftbufs,(5*yd*zd),MPDP,pd.left,tagl,MPI_COMM_WORLD,&V.req[0]);
	MPI_Send_init(V.rightbufs,(5*yd*zd),MPDP,pd.right,tagr,MPI_COMM_WORLD,&V.req[1]);
	MPI_Recv_init(V.leftbufr,(5*yd*zd),MPDP,pd.left,tagr,MPI_COMM_WORLD,&V.req[2]);
	MPI_Recv_init(V.rightbufr,(5*yd*zd),MPDP,pd.right,tagl,MPI_COMM_WORLD,&V.req[3]);
	
	MPI_Send_init(V.dwbufs,(5*xd*zd),MPDP,pd.dw,tagl,MPI_COMM_WORLD,&V.req[4]);
	MPI_Send_init(V.upbufs,(5*xd*zd),MPDP,pd.up,tagr,MPI_COMM_WORLD,&V.req[5]);
	MPI_Recv_init(V.dwbufr,(5*xd*zd),MPDP,pd.dw,tagr,MPI_COMM_WORLD,&V.req[6]);
	MPI_Recv_init(V.upbufr,(5*xd*zd),MPDP,pd.up,tagl,MPI_COMM_WORLD,&V.req[7]);
	
	MPI_Send_init(V.urbufs,9*zd,MPDP,pd.ur,tagl,MPI_COMM_WORLD,&V.req[8]);
	MPI_Recv_init(V.dlbufr,9*zd,MPDP,pd.dl,tagl,MPI_COMM_WORLD,&V.req[9]);
	MPI_Send_init(V.dlbufs,9*zd,MPDP,pd.dl,tagr,MPI_COMM_WORLD,&V.req[10]);
	MPI_Recv_init(V.urbufr,9*zd,MPDP,pd.ur,tagr,MPI_COMM_WORLD,&V.req[11]);

	ch=argv[1][1];
	MPI_Bcast(&ch,1,MPI_CHAR,0,MPI_COMM_WORLD);

	t3=MPI_Wtime();

	switch (ch)
	{
		case 's':
		{
			TM=TSTR;
			V=solid_init(xd,yd,zd,V);
#ifdef RSLIP
			V=slip_ridges_init(xd, yd, zd, V, pd);
#endif
#ifdef PSLIP
			V=slip_posts_init(xd, yd, zd, V, pd);
#endif
			V=init(xd,yd,zd,V);
			ss=0;
			break;
		}
		case 'p':
		{
			if (pd.myrank==0)
				fprintf(stderr,"TRYING TO READ DATA FROM FILES\n");
			cntr=TSTR;

			sprintf(fn,"lbf.%.4d",cntr); 
			MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
			MPI_File_set_view(fh, 0, MPI_DOUBLE, sarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.ftemp, 1, larray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.ftemp %d\n",i);
			
			MPI_File_set_view(fh, (MPI_Offset)(offsett*((MPI_Offset)YDIM)*((MPI_Offset)XDIM)), MPI_DOUBLE, sarray, "native", MPI_INFO_NULL);
			MPI_File_read_all(fh, V.f, 1, larray, &status); 
			MPI_Get_count(&status,MPI_DOUBLE,&i);
			printf("V.f %d\n",i);
			MPI_File_close(&fh);
			
			MPI_Sendrecv(V.ftemp,1,send_to_le,pd.left,4132,V.ftemp,1,recv_from_ri,pd.right,4132,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.f,1,send_to_le,pd.left,4133,V.f,1,recv_from_ri,pd.right,4133,MPI_COMM_WORLD,&status);
			
			MPI_Sendrecv(V.ftemp,1,send_to_dw,pd.dw,4130,V.ftemp,1,recv_from_up,pd.up,4130,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.f,1,send_to_dw,pd.dw,4131,V.f,1,recv_from_up,pd.up,4131,MPI_COMM_WORLD,&status);
			
			MPI_Sendrecv(V.ftemp,1,send_to_dl,pd.dl,4134,V.ftemp,1,recv_from_ur,pd.ur,4134,MPI_COMM_WORLD,&status);
			MPI_Sendrecv(V.f,1,send_to_dl,pd.dl,4135,V.f,1,recv_from_ur,pd.ur,4135,MPI_COMM_WORLD,&status);

			if (pd.myrank==0)
			{
				sprintf(fn,"force.%.4d",cntr);
				sv=fopen(fn,"rb");
				i=fread(&V.force,sizeof(DP),1,sv);
				fclose(sv);
			}
			MPI_Bcast(&V.force,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			printf("V.force=%f myrank=%d was read\n",V.force,pd.myrank);
			
			fprintf(stderr,"end of reading files, myrank=%d\n",pd.myrank);

			V=solid_init(xd,yd,zd,V);
#ifdef RSLIP
			V=slip_ridges_init(xd, yd, zd, V, pd);
#endif
#ifdef PSLIP
			V=slip_posts_init(xd, yd, zd, V, pd);
#endif
			if (pd.myrank==0)
			{
				fprintf(stderr,"END OF READING DATA FROM FILES\n");
				fprintf(stderr,"STARTING FROM %d\n",cntr);
			}
			TM=cntr;
			ss=0;
			break;
		}
		default :
		{
			perror("not a good file name");
			MPI_Finalize();
			return(-1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	V.rho=0.;rhozero=0.;
	for (i=0;i<(xd-1);i++)
	{
		for (j=0;j<(yd-1);j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				rhozero=0.;
				for (a=0;a<19;a++)
					rhozero += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);	
				V.rho+=rhozero;
			}
		}
	}
	V.rho /= ((double)((xd-1)*(yd-1)*XDIM*YDIM*(zd-2)));
	MPI_Allreduce(&V.rho,&rhozero,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    	V.rho=rhozero;
	printf("Average Density in the channel at start: %f  myrank=%d\n",V.rho,pd.myrank);

	t4=MPI_Wtime();

 	for (ss;ss<(up+1);ss++)
	{
		V=encalc(xd, yd , zd, V, Gx,TM,ss,ubulk, cart_grid);
		V=wshearstr(xd, yd, zd, V, Gx, tau, TM, ss, cart_grid);
		if (!(ss%(int)utime))
		{
			if (!(ss%(10*(int)utime)))
			{
				sprintf(fn,"lbf.%.4ld",(TM+(ss/(int)utime)));
				MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
				MPI_File_set_view(fh, 0, MPI_DOUBLE, sarray, "native", info);
				MPI_File_write_all(fh, V.ftemp, 1, larray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.ftemp %d\n",i);
				MPI_File_set_view(fh, (MPI_Offset)(offsett*((MPI_Offset)YDIM)*((MPI_Offset)XDIM)), MPI_DOUBLE, sarray, "native", info);
				MPI_File_write_all(fh, V.f, 1, larray, &status); 
				MPI_Get_count(&status,MPI_DOUBLE,&i);
				printf("V.f %d\n",i);
				MPI_File_close(&fh);
				
				if (pd.myrank==0)
				{
					sprintf(fn,"force.%.4ld",(TM+(ss/(int)utime)));
					sv=fopen(fn,"wb");
					i=fwrite(&V.force,sizeof(DP),1,sv);
					fclose(sv);
				}
				MPI_Bcast(&i,1,MPI_INT,0,MPI_COMM_WORLD);
				printf("writing V.force=%f num written=%d myrank=%d   V.rho=%f\n",V.force,i,pd.myrank,V.rho);
			}
			V=vel_print(xd,yd,zd,V,(TM+(ss/(int)utime)),cart_grid);
			V=statistics(xd,yd,zd,V,pd,(TM+(ss/(int)utime)));
		}
		V=xcomputations(xd,yd,zd,V,pd,Gx,Gy,Gz,rhozero,tau,ubulk);
	}

	t2=MPI_Wtime();
	if (!pd.myrank)
	{
		tm=fopen("report.code","w");
		V=output_twoD(xd,yd,zd,V,C,dx,Fx,rhozero,tau);
		fprintf(tm,"xd=%d,yd=%d,zd+2=%d,Reb=%f,tau=%f,G=%.10f,Ub=%f,CFLnumber=%f\n",xd,yd,zd,Reb,tau,Fx,ubulk,CFL);
		fprintf(tm,"number of processors: %d\n",pd.numproc);
		fprintf(tm,"Executed File's Name: %s\n",__FILE__);
		fprintf(tm,"Compilation Option: xcut\n");
#ifdef RSLIP
		fprintf(tm,"slip option: RSLIP\n");
#endif
#ifdef PSLIP
		fprintf(tm,"slip option: PSLIP\n");
#endif
		fprintf(tm,"Total computational time is: %f s\n",t2-t4);
		fprintf(tm,"Unit non-dimensional time is: %f iterations\n",utime);
		fprintf(tm,"Starting time: %d\n",TSTR);
		fprintf(tm,"Current Time is: %d\n",TEND);
		fclose(tm);
	}
	
	free(V.f);
	free(V.ftemp);
	free(V.leftbufr);
	free(V.leftbufs);
	free(V.rightbufr);
	free(V.rightbufs);
	free(V.dwbufr);
	free(V.dwbufs);
	free(V.upbufr);
	free(V.upbufs);
	free(V.urbufr);
	free(V.dlbufs);
	
	free(V.s);

	for (i=0;i<12;i++)
	  MPI_Request_free(&V.req[i]);
	
	MPI_Type_free(&block);
	MPI_Type_free(&larray);
	MPI_Type_free(&sarray);
	MPI_Type_free(&send_to_dw);
	MPI_Type_free(&recv_from_up);
	MPI_Type_free(&send_to_le);
	MPI_Type_free(&recv_from_ri);
	MPI_Type_free(&send_to_dl);
	MPI_Type_free(&recv_from_ur);
	
	MPI_Info_free(&info);

	MPI_Finalize();

	return(0);
}

