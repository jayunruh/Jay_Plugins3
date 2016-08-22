/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jalgs.*;
import jalgs.jfit.*;
import ij.text.*;
import jguis.*;

public class fit_traj_gaussian_gridsearch_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	int xpts,ypts;
	float[] data;
	gausfunc gf;

	public void run(String arg) {
		ImagePlus[] imp=jutils.selectImages(false,2,new String[]{"Image","Trajectory"});
		if(imp==null) return;
		int width=imp[0].getWidth();
		int height=imp[0].getHeight();
		float[][] charrays=new float[2][];
		ImageWindow iw=imp[1].getWindow();
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
		if(sel<0) sel=0;
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		String[] annotations=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");
		int frameoffset=0;
		if(annotations!=null) frameoffset=(int)Float.parseFloat(annotations[sel]);
		gf=new gausfunc();
		gridfit fitclass=new gridfit(this);
		fitclass.output=false;
		xpts=10; ypts=10;
		boolean findmax=true;
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Find_Z_Max",findmax);
		gd2.addNumericField("X_Size (pixels)",xpts,0);
		gd2.addNumericField("Y_Size (pixels)",ypts,0);
		gd2.addNumericField("xc_min (rel)",-2.0,5,10,null);
		gd2.addNumericField("xc_max (rel)",2.0,5,10,null);
		gd2.addNumericField("yc_min (rel)",-2.0,5,10,null);
		gd2.addNumericField("yc_max (rel)",2.0,5,10,null);
		gd2.addNumericField("stdev_min",0.5,5,10,null);
		gd2.addNumericField("stdev_max",(double)xpts*0.25,5,10,null);
		gd2.addCheckbox("Calibrate",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		findmax=gd2.getNextBoolean();
		xpts=(int)gd2.getNextNumber();
		ypts=(int)gd2.getNextNumber();
		double xcmin=gd2.getNextNumber();
		double xcmax=gd2.getNextNumber();
		double ycmin=gd2.getNextNumber();
		double ycmax=gd2.getNextNumber();
		double stdevmin=gd2.getNextNumber();
		double stdevmax=gd2.getNextNumber();
		boolean calibrate=gd2.getNextBoolean();

		int currz=imp[0].getZ();
		int currt=imp[0].getT();
		int currc=imp[0].getC();
		int slices=imp[0].getNSlices();
		int frames=imp[0].getNFrames();
		int nchannels=imp[0].getNChannels();
		int maxslice=currz-1;
		float[][] fits=new float[npts[sel]][];
		float[][] datas=new float[npts[sel]][];
		//need to loop over the trajectory frames here
		for(int i=0;i<npts[sel];i++){
			int centerx=(int)(xvals[sel][i]+0.5f);
			int centery=(int)(yvals[sel][i]+0.5f);
			float xoff=(float)(centerx-xpts/2);
			float yoff=(float)(centery-ypts/2);
			Rectangle r=new Rectangle(centerx-xpts/2,centery-ypts/2,xpts,ypts);
			Object temp=jutils.get3DSlice(imp[0].getStack(),frameoffset+i,maxslice,currc-1,frames,slices,nchannels);
			if(findmax){
				Object[] temp1=jutils.get3DZSeries(imp[0].getStack(),currc-1,frameoffset+i,frames,slices,nchannels);
				float maxint=jstatistics.getstatistic("Avg",temp1[0],width,height,r,null);
				maxslice=0;
				for(int j=1;j<slices;j++){
					float tempint=jstatistics.getstatistic("Avg",temp1[j],width,height,r,null);
					if(tempint>maxint){
						maxint=tempint;
						maxslice=j;
					}
				}
				temp=temp1[maxslice];
			}
			data=algutils.convert_arr_float(algutils.get_region(temp,centerx,centery,xpts,ypts,width,height));
			//datas[i]=data.clone();
			//float[] xymax=getxymax(data,xpts,ypts);
			double[] params=new double[3];
			double[][] constraints=new double[3][3];
			constraints[0][0]=stdevmin; constraints[1][0]=stdevmax; constraints[2][0]=0.1;
			constraints[0][1]=centerx-xoff+xcmin; constraints[1][1]=centerx-xoff+xcmax; constraints[2][1]=0.25;
			constraints[0][2]=centery-yoff+ycmin; constraints[1][2]=centery-yoff+ycmax; constraints[2][2]=0.25;
			double[] stats=new double[2];
			int[] fixes={0,0,0};
			fits[i]=fitclass.fitdata(params,fixes,constraints,data,null,stats);
			double[] ampoffset=(new linleastsquares()).get_amp_offset(get_gaussian(params),data,true);
			if(calibrate){
				double psize=jutils.get_psize(imp[0]);
				params[0]*=psize;
				params[1]*=psize;
				params[2]*=psize;
				xoff*=(float)psize;
				yoff*=(float)psize;
			}
			TextWindow tw=jutils.selectTable("Gaussian_Output");
			if(tw==null) tw=new TextWindow("Gaussian_Output","title\tbaseline\tamplitude\tstdev\txc\tyc\tzslice\tframe\tsel_traj\tc2","",400,200);
			tw.append(imp[0].getTitle()+"\t"+(float)ampoffset[1]+"\t"+(float)ampoffset[0]+"\t"+(float)params[0]+"\t"+((float)params[1]+xoff)+"\t"+((float)params[2]+yoff)+"\t"+maxslice+"\t"+(frameoffset+i+1)+"\t"+sel+"\t"+stats[1]+"\n");
			IJ.showProgress(i,npts[sel]);
		}
		new ImagePlus("fits",jutils.array2stack(fits,xpts,ypts)).show();
		//new ImagePlus("datas",jutils.array2stack(datas,xpts,ypts)).show();
	}

	public float[] getxymax(float[] image,int width,int height){
		float max=image[0];
		int xmax=0;
		int ymax=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]>max){
					max=image[j+i*width];
					xmax=j;
					ymax=i;
				}
			}
		}
		return new float[]{max,(float)xmax,(float)ymax};
	}

	public double[] fitfunc(double[] params)
	{
		//the params list is stdev,xc,yc
		double[] func=get_gaussian(params);
		double[] ampoffset=(new linleastsquares()).get_amp_offset(func,data,true);
		for(int i=0;i<(xpts*ypts);i++){func[i]*=ampoffset[0]; func[i]+=ampoffset[1];}
		return func;
	}

	public double[] get_gaussian(double[] params){
		float[] func= gf.get_2D_func2(-params[1],xpts,1.0,-params[2],ypts,1.0,params[0]);
		return algutils.convert_arr_double(func);
		/*double[] func=new double[xpts*ypts];
		for(int i=0;i<ypts;i++){
			double ydiff=((double)i-params[2]);
			for(int j=0;j<xpts;j++){
				double xdiff=((double)j-params[1]);
				func[j+i*xpts]= Math.exp(((-0.5)*(xdiff*xdiff))/(params[0]*params[0]));
				func[j+i*xpts]*= Math.exp(((-0.5)*(ydiff*ydiff))/(params[0]*params[0]));
			}
		}
		return func;*/
	}

	public void showresults(String results){
		IJ.log(results);
	}


}
