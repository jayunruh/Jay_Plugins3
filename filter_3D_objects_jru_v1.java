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
import ij.plugin.frame.*;
import jalgs.*;
import jalgs.jseg.*;

public class filter_3D_objects_jru_v1 implements PlugIn, gui_interface {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean clearedges=true;
		gd.addCheckbox("Clear_Edges?",clearedges);
		int minarea=10;
		gd.addNumericField("Min Size (voxels)",minarea,0);
		int maxarea=100;
		gd.addNumericField("Max Size (voxels)",maxarea,0);
		boolean newimage=true;
		gd.addCheckbox("Create_New_Image?",newimage);
		boolean outindexed=false;
		gd.addCheckbox("Output_Indexed_Image?",outindexed);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		clearedges=gd.getNextBoolean();
		minarea=(int)gd.getNextNumber();
		maxarea=(int)gd.getNextNumber();
		newimage=gd.getNextBoolean();
		outindexed=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		if(slices==1){slices=frames; frames=1;}
		findblobs3D fb=new findblobs3D(width,height,slices,this);
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<frames;i++){
			byte[][] data=new byte[slices][];
			for(int j=0;j<slices;j++) data[j]=(byte[])stack.getPixels(i*slices+j+1);
			float[][] objects=fb.dofindblobs(data);
			if(clearedges) fb.clear_edges(objects,false);
			int[] filter={minarea,maxarea};
			fb.filter_area(objects,filter,true);
			if(newimage || outindexed){
				if(outindexed) for(int j=0;j<slices;j++) retstack.addSlice("",objects[j]);
				else{
					byte[][] temp=fb.tobinary(objects,true);
					for(int j=0;j<slices;j++) retstack.addSlice("",temp[j]);
				}
			} else {
				byte[][] temp=fb.tobinary(objects,true);
				for(int j=0;j<slices;j++) stack.setPixels(temp[j],i*slices+j+1);
				//data=fb.tobinary(objects,true);
			}
		}
		if(newimage || outindexed){
			new ImagePlus(imp.getTitle()+" filtered",retstack).show();
		} else {
			imp.setStack(null,stack);
			imp.updateAndDraw();
		}
	}

	public void showMessage(String message){}

	public void showProgress(int currpos,int finalpos){IJ.showProgress(currpos,finalpos);}

}
