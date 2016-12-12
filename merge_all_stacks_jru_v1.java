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
import jguis.*;

public class merge_all_stacks_jru_v1 implements PlugIn {

	public void run(String arg) {
		//assume all images have same number of channels and slices
		int[] wList=WindowManager.getIDList();
		ImagePlus imp=WindowManager.getImage(wList[0]);
		ImageStack stack=imp.getStack();
		int channels=imp.getNChannels();
		int slices=imp.getNSlices();
		int totsize=0;
		int[] widths=new int[wList.length];
		int[] heights=new int[wList.length];
		//int[] channels=new int[wList.length];
		//int[] slices=new int[wList.length];
		int maxwidth=0;
		int maxheight=0;
		float psize=(float)jutils.get_psize(imp);
		//int maxchans=0;
		//int maxslices=0;
		for(int i=0;i<wList.length;i++){
			ImagePlus imp2=WindowManager.getImage(wList[i]);
			if(imp2!=null){
				widths[i]=imp2.getWidth();
				heights[i]=imp2.getHeight();
				//channels[i]=imp2.getNChannels();
				//slices[i]=imp2.getNSlices();
				if(widths[i]>maxwidth) maxwidth=widths[i];
				if(heights[i]>maxheight) maxheight=heights[i];
				//if(channels[i]>maxchans) maxchans=channels[i];
				//if(slices[i]>maxwidth) maxslices=slices[i];
			}
		}
		ImageStack newstack=new ImageStack(maxwidth,maxheight);
		for(int i=0;i<wList.length;i++){
			ImagePlus imp2=WindowManager.getImage(wList[i]);
			if(imp2!=null){
				ImageStack stack2=imp2.getStack();
				int size=stack2.getSize();
				totsize+=size;
				for(int j=0;j<size;j++){
					Object pixels=stack2.getPixels(1);
					//float[] copied=new float[maxwidth*maxheight];
					int type=algutils.get_array_type(pixels);
					Object copied=algutils.create_array(maxwidth*maxheight,type);
					float xshift=0.5f*(float)(maxwidth-widths[i]);
					float yshift=0.5f*(float)(maxheight-heights[i]);
					interpolation.shift_copy_image(pixels,widths[i],heights[i],copied,maxwidth,maxheight,xshift,yshift);
					newstack.addSlice(imp2.getTitle(),copied);
					stack2.deleteSlice(1);
				}
				imp2.changes=false;
				imp2.close();
			}
		}
		int framesize=channels*slices;
		int nframes=(int)(totsize/framesize);
		ImagePlus imp3=new ImagePlus("Merged Stacks",newstack);
		imp3.setOpenAsHyperStack(true);
		imp3.setDimensions(channels,slices,nframes);
		jutils.set_psize(imp3,(double)psize);
		if(channels>1) new CompositeImage(imp3,CompositeImage.COLOR).show();
		else imp3.show();
	}

}
