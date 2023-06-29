#include "definitions.h"

POINTER vel_print(int xd, int yd, int zd, POINTER V, long cntr, MPI_Comm cart_grid)
{
	int i,j,k,a,q;
	char fl[40];
	double ux,uy,uz,dens,*tmp,*tmpd;
	MPI_File fh,fhd;
	MPI_Offset offset,offsetd;
	MPI_Status status,statusd;
	
	int coord[2], rank, size;
	MPI_Datatype block,cblock,sarray,d_sarray;
	MPI_Datatype cdblock;
	MPI_Info info;
	int gsizes[2],lsizes[2],psizes[2],array_start[2],memsizes[2];

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Cart_coords(cart_grid, rank, 2, coord);
	
	MPI_Info_create(&info);
	MPI_Info_set(info, "striping_factor", STRIPE_COUNT);
	MPI_Info_set(info, "cb_nodes", CB_NODES_FCT);
//	MPI_Info_set(info, "striping_unit", STRIPE_SIZE);
	MPI_Info_set(info, "romio_cb_write", "enable");
//	MPI_Info_set(info, "romio_ds_write", "disable");
	
	MPI_Type_contiguous(3,MPI_DOUBLE,&block);
	MPI_Type_contiguous(zd,MPI_DOUBLE,&cdblock);
	MPI_Type_contiguous(zd,block,&cblock);
	
	MPI_Type_commit(&cdblock);
	MPI_Type_commit(&cblock);
	
	gsizes[0] = (xd-1)*XDIM;    gsizes[1] = (yd-1)*YDIM;     
	psizes[0] = XDIM;    psizes[1] = YDIM;
	lsizes[0] = xd-1; lsizes[1] = yd-1;
	array_start[0] = coord[0] * lsizes[0];
	array_start[1] = coord[1] * lsizes[1];
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, cblock, &sarray);
	MPI_Type_commit(&sarray);
	MPI_Type_create_subarray(2, gsizes, lsizes, array_start, MPI_ORDER_C, cdblock, &d_sarray);
	MPI_Type_commit(&d_sarray);
	
	
	tmp=(double *)calloc(3*zd*(yd-1)*(xd-1),sizeof(double));
	tmpd=(double *)calloc(zd*(yd-1)*(xd-1),sizeof(double));// it should be xdc-1
	
	if (tmp==NULL)
		fprintf(stderr,"Unable to allocte memory for vel_print function tmp\n error code: %s\n",strerror(errno));

	if (tmpd==NULL)
		fprintf(stderr,"Unable to allocte memory for vel_print function tmpd\n error code: %s\n",strerror(errno));

	sprintf(fl,"vel.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fh);
	MPI_File_set_view(fh, 0, MPI_DOUBLE, sarray, "native", info);

	sprintf(fl,"den.%.4d",((int)cntr));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,info,&fhd);
	MPI_File_set_view(fhd, 0, MPI_DOUBLE, d_sarray, "native", info);


	for (i=0;i<(xd-1);i++)
	{
		for (j=0;j<(yd-1);j++)
		{
			for (k=0;k<zd;k++)
			{
				ux=0.; uy=0.; uz=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					ux += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					uy += E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					uz += E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					dens += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				}
				if (dens>0.)
				{
					ux /= dens;
					uy /= dens;
					uz /= dens;
				}
				else
				{
					ux =0.;
					uy =0.;
					uz =0.;
					dens=0.;
				}

				*(tmp+i*(yd-1)*zd*3+j*zd*3+k*3+0)=ux;
				*(tmp+i*(yd-1)*zd*3+j*zd*3+k*3+1)=uy;
				*(tmp+i*(yd-1)*zd*3+j*zd*3+k*3+2)=uz;
 				*(tmpd+i*(yd-1)*zd+j*zd+k)=dens;
			}
		}
	}
 	MPI_File_write_all(fh,tmp,(xd-1)*(yd-1)*zd*3,MPI_DOUBLE,&status);
	MPI_File_write_all(fhd,tmpd,(xd-1)*(yd-1)*zd,MPI_DOUBLE,&statusd);
	
	free(tmp);
	free(tmpd);
	
	MPI_File_close(&fh);
	MPI_File_close(&fhd);

	MPI_Type_free(&cdblock);
	MPI_Type_free(&cblock);
	
	MPI_Type_free(&sarray);
	MPI_Type_free(&d_sarray);
	
	MPI_Info_free(&info);
		
	return V;
}