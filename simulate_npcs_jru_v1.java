/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import jalgs.jfit.*;
import jalgs.jsim.*;

public class simulate_npcs_jru_v1 implements PlugIn {
	//this plugin simulates random points on the surface of a sphere convolved with a 3D psf
	//the image is 128 x 128 x 128

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Sphere_Radius (nm)",1000.0,5,15,null);
		gd.addNumericField("Number_of_points",50,0);
		gd.addNumericField("PSF_FWHM (nm)",200.0,5,15,null);
		gd.addNumericField("PSF_Z_FWHM (nm)",600.0,5,15,null);
		gd.addNumericField("Pixel_Size (nm)",40.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float srad=(float)gd.getNextNumber();
		int npts=(int)gd.getNextNumber();
		float xystdev=((float)gd.getNextNumber())/2.35f;
		float zstdev=((float)gd.getNextNumber())/2.35f;
		float psize=(float)gd.getNextNumber();
		float[][] spherepts=makeSphere(npts,srad,64.0f*psize,64.0f*psize,64.0f*psize);
		plotPoints(spherepts);
		Object[] stack=new Object[128];
		for(int i=0;i<128;i++){
			stack[i]=new float[128*128];
		}
		drawPoints(stack,spherepts,psize,psize,xystdev,zstdev,128,128);
		ImagePlus imp=new ImagePlus("Sim NPCs",jutils.array2stack(stack,128,128));
		jutils.set_psize(imp,psize/1000.0f);
		imp.show();
	}

	public void plotPoints(float[][] coords){
		float[][] xpts=new float[coords.length][1];
		float[][] ypts=new float[coords.length][1];
		float[][] zpts=new float[coords.length][1];
		for(int i=0;i<coords.length;i++){
			xpts[i][0]=coords[i][0]/1000.0f;
			ypts[i][0]=coords[i][1]/1000.0f;
			zpts[i][0]=coords[i][2]/1000.0f;
		}
		Traj3D traj=new Traj3D("x","y","z",xpts,ypts,zpts,null);
		int[] shapes=traj.getShapes();
		for(int i=0;i<coords.length;i++){
			shapes[i]=1;
		}
		new PlotWindow3D("Sphere Positions",traj).draw();
	}

	public float[][] makeSphere(int nparticles,float radius,float xc, float yc, float zc){
		float[][] coords=new float[nparticles][3];
		rngs random=new rngs();
		for(int i=0;i<nparticles;i++){
			double[] dcoords=random.random_sphere(radius);
			coords[i][0]=(float)dcoords[0]+xc;
			coords[i][1]=(float)dcoords[1]+yc;
			coords[i][2]=(float)dcoords[2]+zc;
		}
		return coords;
	}

	public void drawPoints(Object[] stack,float[][] coords,float psize,float zsize,float xystdev,float zstdev,int width,int height){
		gausfunc gf=new gausfunc();
		float xystdevpix=xystdev/psize;
		float zstdevpix=zstdev/zsize;
		for(int i=0;i<coords.length;i++){
			float xpos=coords[i][0]/psize;
			float ypos=coords[i][1]/psize;
			float zpos=coords[i][2]/zsize;
			gf.draw_3D_func(stack,xpos,ypos,zpos,width,height,xystdevpix,zstdevpix,1.0f);
		}
		return;
	}

}
