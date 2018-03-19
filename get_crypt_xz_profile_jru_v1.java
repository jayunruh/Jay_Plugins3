/*******************************************************************************
 * Copyright (c) 2017 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class get_crypt_xz_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		float zratio=(float)jutils.get_zratio(imp);
		float zoff=50f;
		int zsize=250;
		float expansion=1.5f;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addNumericField("Z_Offset (pixel units)",zoff,5,15,null);
		gd.addNumericField("Z_Size (pixel units)",zsize,0);
		gd.addNumericField("XY_Size (fraction of diameter)",expansion,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		zratio=(float)gd.getNextNumber();
		zoff=(float)gd.getNextNumber();
		zsize=(int)gd.getNextNumber();
		expansion=(float)gd.getNextNumber();
		ImageStack stack=imp.getStack();
		int nchans=imp.getNChannels();
		int nslices=imp.getNSlices();
		Object[][] stack2=new Object[nchans][nslices];
		Object[] temp=jutils.stack2array(stack);
		for(int i=0;i<nslices;i++){
			for(int j=0;j<nchans;j++){
				stack2[j][i]=temp[j+i*nchans];
			}
		}
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) {
			IJ.error("need crypt neck ellipses and vertex points in roi manager");
			return;
		}
		//each pair of items in the roi manager will be an ellipse around the crypt neck and then a point at the crypt vertex
		//need to transform each crypt so it is oriented vertically
		//then crop it (use max ellipse diameter) and create the radial profile in all channels
		Roi[] rois=rman.getRoisAsArray();
		int ncrypts=rois.length/2;
		for(int i=0;i<ncrypts;i++){
			IJ.log("analyzing crypt "+(i+1));
			//for each crypt get the centroid of the ellipse and the vertex and the longest dimension of the ellipse
			Roi vertexroi=rois[2*i+1];
			Rectangle r=vertexroi.getBounds();
			Roi neckroi=rois[2*i];
			float[] vertex={r.x,r.y,zratio*(float)(vertexroi.getZPosition()-1)};
			IJ.log("vertex pos = \t"+table_tools.print_float_array(vertex));
			double[] params=((EllipseRoi)neckroi).getParams();
			float[] neck={0.5f*(float)(params[0]+params[2]),0.5f*(float)(params[1]+params[3]),zratio*(float)(neckroi.getZPosition()-1)};
			IJ.log("neck center = \t"+table_tools.print_float_array(neck));
			float maxd=(float)Math.sqrt((params[2]-params[0])*(params[2]-params[0])+(params[3]-params[1])*(params[3]-params[1]));
			//for the transformation, need to rotate about the cross-product between the crypt vector and the z axis
			//rotation angle is given by the dot product between the crypt vector and the z axis
			float[] cryptvec={(neck[0]-vertex[0]),(neck[1]-vertex[1]),(neck[2]-vertex[2])};
			cryptvec=measure_object.norm_vector(cryptvec);
			float[] zvec={0.0f,0.0f,1.0f};
			float angle=measure_object.get_inner_angle(zvec,cryptvec); //order?
			IJ.log("angle = \t"+angle);
			angle=0.0f;
			//now get the cross product
			float[] crossprod=measure_object.crossProd(zvec,cryptvec);
			crossprod=measure_object.norm_vector(crossprod);
			IJ.log("rot vector = \t"+table_tools.print_float_array(crossprod));
			float[][] rotmat=measure_object.getRotationMatrix(crossprod,angle);
			int rotsize=(int)(expansion*maxd);
			int rsize=(int)(0.5f*(float)rotsize-0.5f);
			float[][] xzstack=new float[nchans][];
			for(int j=0;j<nchans;j++){
				Object[] chanstack=stack2[j];
				float[][] rotated=getCryptImage(stack2[j],width,height,vertex,zratio,rotmat,maxd,zoff,zsize,expansion);
				//new ImagePlus("rotated",jutils.array2stack(rotated,rotsize,rotsize)).show();
				xzstack[j]=getxzProfile(rotated,rotsize,rotsize,rsize);
			}
			jutils.create_hyperstack("crypt "+i+" xzprofile",jutils.array2stack(xzstack,rsize*2-1,zsize),imp,1,1,nchans).show();
		}
	}

	public float[] getxzProfile(float[][] stack,int width,int height,int rsize){
		int newsize=2*rsize-1;
		float[] xzprof=new float[newsize*stack.length];
		for(int i=0;i<stack.length;i++){
			float[] circavg=interpolation.circavg(stack[i],width,height,rsize,width/2,height/2,true);
			System.arraycopy(circavg,0,xzprof,i*newsize,newsize);
		}
		return xzprof;
	}

	public float[][] getCryptImage(Object[] stack,int width,int height,float[] vertex,float zratio,float[][] rotmat,float maxd,float zshift,int zsize,float expansion){
		//now build a rotated crypt by rotating each voxel to its position in the original image
		int rotsize=(int)(expansion*maxd);
		//int rotheight=4*rotsize;
		int rotheight=zsize;
		float[][] rotated=new float[rotheight][rotsize*rotsize];
		for(int i=0;i<rotheight;i++){
			float zpos=(float)(i-zshift);
			for(int j=0;j<rotsize;j++){
				float ypos=(float)(j-rotsize/2);
				for(int k=0;k<rotsize;k++){
					float xpos=(float)(k-rotsize/2);
					//if(i==0 && j==0 && k==0) IJ.log(""+xpos+" , "+ypos+" , "+zpos);
					//multiply by the rotation matrix to tranform int the old coordinates
					float xoff=rotmat[0][0]*xpos+rotmat[0][1]*ypos+rotmat[0][2]*zpos;
					float yoff=rotmat[1][0]*xpos+rotmat[1][1]*ypos+rotmat[1][2]*zpos;
					float zoff=rotmat[2][0]*xpos+rotmat[2][1]*ypos+rotmat[2][2]*zpos;
					//correct for anisotropic resolution
					zoff/=zratio;
					//add the vertex position back
					xoff+=vertex[0]; yoff+=vertex[1]; zoff+=(vertex[2]/zratio);
					//if(i==0 && j==0 && k==0) IJ.log(""+xoff+" , "+yoff+" , "+zoff);
					//finally interpolate the original image at these points
					rotated[i][k+j*rotsize]=interpolation.interp3D(stack,width,height,xoff,yoff,zoff);
				}
			}
		}
		return rotated;
	}

}
