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
import jguis.*;
import jalgs.*;

public class measure_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we measure trajectories on a stack
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Image","Trajectory"});
		if(imps==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Statistic",jstatistics.stats,"Avg");
		gd.addNumericField("Measure_Radius",5.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		float rad=(float)gd.getNextNumber();
		ImageWindow iw=imps[1].getWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		String[] starts=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");

		ImageStack stack=imps[0].getStack();
		int width=imps[0].getWidth(); int height=imps[0].getHeight();
		int frames=stack.getSize(); //assume single channel for now
		//IJ.log(""+npts.length);
		int totlength=0;
		for(int i=0;i<npts.length;i++) totlength+=npts[i];
		float[][] stats=new float[totlength][6];
		int counter=0;
		//stats are id,start,frame,x,y,stat
		for(int i=0;i<npts.length;i++){
			int start=0;
			if(starts!=null){
				start=(int)Float.parseFloat(starts[i]);
			}
			for(int j=start;j<(start+npts[i]);j++){
				float[] frame=algutils.convert_arr_float2(stack.getPixels(j+1));
				float[] circ=getCircleVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],rad);
				stats[counter][0]=i;
				stats[counter][1]=start;
				stats[counter][2]=j-start;
				stats[counter][3]=xvals[i][j-start];
				stats[counter][4]=yvals[i][j-start];
				stats[counter][5]=jstatistics.getstatistic(stat,circ,null);
				counter++;
			}
		}
		//now output the stats
		table_tools.create_table("Traj Stats",stats,new String[]{"id","start","frame","x","y",stat});
	}

	public float[] getCircleVals(float[] image,int width,int height,float xc,float yc,float rad){
		int size=(int)(2.0f*rad+1.0f);
		float[] square=new float[size*size];
		int startx=(int)(xc-rad);
		int starty=(int)(yc-rad);
		int counter=0;
		for(int i=starty;i<(starty+size);i++){
			float y=yc-(float)i;
			for(int j=startx;j<(startx+size);j++){
				float x=xc-(float)j;
				float dist2=x*x+y*y;
				if(dist2<=rad && i>=0 && i<height && j>=0 && j<width){
					square[counter]=image[j+i*width];
					counter++;
				}
			}
		}
		return (float[])algutils.get_subarray(square,0,counter);		
	}

	public PolygonRoi traj2roi(float[] xvals,float[] yvals,int npts){
		int[] xvals2=new int[npts];
		int[] yvals2=new int[npts];
		for(int i=0;i<npts;i++){
			xvals2[i]=(int)xvals[i];
			yvals2[i]=(int)yvals[i];
		}
		return new PolygonRoi(xvals2,yvals2,npts,Roi.POLYLINE);
	}

	public int[][] traj2int(float[] xvals,float[] yvals,int npts){
		int[] xvals2=new int[npts];
		int[] yvals2=new int[npts];
		for(int i=0;i<npts;i++){
			xvals2[i]=(int)xvals[i];
			yvals2[i]=(int)yvals[i];
		}
		return new int[][]{xvals2,yvals2};
	}

}
