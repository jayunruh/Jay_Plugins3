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
		gd.addNumericField("Measure_Radius_Z (optional)",5.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		float rad=(float)gd.getNextNumber();
		float zrad=(float)gd.getNextNumber();
		ImageWindow iw=imps[1].getWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		String[] starts=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");
		boolean threed=jutils.is3DPlot(iw);
		int[] npts=null;
		if(threed) npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
		else npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[][] zvals=null;
		if(threed) zvals=((float[][][])jutils.runPW4VoidMethod(iw,"getZValues"))[0];

		ImageStack stack=imps[0].getStack();
		int width=imps[0].getWidth(); int height=imps[0].getHeight();
		int channels=imps[0].getNChannels();
		int currchan=imps[0].getC();
		int slices=imps[0].getNSlices();
		int frames=imps[0].getNFrames();
		if(threed && slices==1){slices=frames; frames=1;}
		//int frames=stack.getSize(); //assume single channel for now
		//IJ.log(""+npts.length);
		int totlength=0;
		for(int i=0;i<npts.length;i++) totlength+=npts[i];
		float[][] stats=new float[totlength][6];
		if(threed) stats=new float[totlength][7];
		int counter=0;
		//stats are id,start,frame,x,y,(z),stat
		for(int i=0;i<npts.length;i++){
			int start=0;
			if(starts!=null){
				start=(int)Float.parseFloat(starts[i]);
			}
			for(int j=start;j<(start+npts[i]);j++){
				float[] circ=null;
				if(threed){
					Object[] frame=jutils.get3DZSeries(stack,currchan-1,j,frames,slices,channels);
					circ=getSphereVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],zvals[i][j-start],rad,zrad);
				} else { 
					Object frame=jutils.get3DSlice(stack,j,0,currchan-1,frames,slices,channels);
					circ=getCircleVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],rad);
				}
				stats[counter][0]=i;
				stats[counter][1]=start;
				stats[counter][2]=j-start;
				stats[counter][3]=xvals[i][j-start];
				stats[counter][4]=yvals[i][j-start];
				if(threed){
					stats[counter][5]=zvals[i][j-start];
					stats[counter][6]=jstatistics.getstatistic(stat,circ,null);
				} else {
					stats[counter][5]=jstatistics.getstatistic(stat,circ,null);
				}
				counter++;
			}
		}
		//now output the stats
		if(threed) table_tools.create_table("Traj Stats",stats,new String[]{"id","start","frame","x","y","z",stat});
		else table_tools.create_table("Traj Stats",stats,new String[]{"id","start","frame","x","y",stat});
	}

	public float[] getCircleVals(Object image,int width,int height,float xc,float yc,float rad){
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
					int pos=j+i*width;
					Object temp=image;
					if(temp instanceof float[]) square[counter]=((float[])temp)[pos];
					if(temp instanceof short[]) square[counter]=(float)(((short[])temp)[pos]&0xffff);
					if(temp instanceof byte[]) square[counter]=(float)(((byte[])temp)[pos]&0xff);
					counter++;
				}
			}
		}
		return (float[])algutils.get_subarray(square,0,counter);		
	}

	public float[] getSphereVals(Object[] image,int width,int height,float xc,float yc,float zc,float rad,float zrad){
		int size=(int)(2.0f*rad+1.0f);
		int zsize=(int)(2.0f*zrad+1.0f);
		int slices=image.length;
		float[] square=new float[size*size*zsize];
		int startx=(int)(xc-rad);
		int starty=(int)(yc-rad);
		int startz=(int)(zc-zrad);
		float zratio=zrad/rad;
		int counter=0;
		for(int i=starty;i<(starty+size);i++){
			float y=yc-(float)i;
			for(int j=startx;j<(startx+size);j++){
				float x=xc-(float)j;
				for(int k=startz;k<(startz+zsize);k++){
					float z=zc-(float)k;
					float dist2=x*x+y*y+z*z/(zratio*zratio);
					if(dist2<=rad && i>=0 && i<height && j>=0 && j<width && k>=0 && k<slices){
						Object temp=image[k];
						int pos=j+i*width;
						if(temp instanceof float[]) square[counter]=((float[])temp)[pos];
						if(temp instanceof short[]) square[counter]=(float)(((short[])temp)[pos]&0xffff);
						if(temp instanceof byte[]) square[counter]=(float)(((byte[])temp)[pos]&0xff);
						counter++;
					}
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
