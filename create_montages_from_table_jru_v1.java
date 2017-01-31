/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import ij.text.*;
import java.util.*;
import jguis.*;
import ij.io.*;
import jalgs.jseg.*;

public class create_montages_from_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addNumericField("Number_Of_Columns",2,0);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		int ncols=(int)gd2.getNextNumber();
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<ncols;i++){
			gd.addChoice("Image"+(i+1)+"_Column",col_labels,col_labels[0]);
		}
		gd.addChoice("Slice_Name_Column",col_labels,col_labels[0]);
		gd.addNumericField("Bin_By",2,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int[] selcols=new int[ncols];
		for(int i=0;i<ncols;i++){
			selcols[i]=gd.getNextChoiceIndex();
		}
		int scol=gd.getNextChoiceIndex();
		int binby=(int)gd.getNextNumber();
		DirectoryChooser dc=new DirectoryChooser("Choose Directory");
		String dir=dc.getDirectory();
		if(dir==null) return;
		ImageStack stack=null;
		int nchans=1;
		int nslices=1;
		int width=1;
		int height=1;
		for(int i=0;i<listtable.size();i++){
			ImagePlus[] imps=new ImagePlus[ncols];
			for(int j=0;j<ncols;j++){
				imps[j]=IJ.openImage(dir+listtable.get(i).get(selcols[j]));
				if(imps[j]==null) break;
			}
			if(imps[ncols-1]==null) continue;
			String name=listtable.get(i).get(scol);
			ImageStack[] stacks=new ImageStack[ncols];
			for(int j=0;j<ncols;j++) stacks[j]=imps[j].getStack();
			if(i==0){
				width=imps[0].getWidth();
				height=imps[0].getHeight();
				stack=new ImageStack(2*width/binby,height/binby);
				nchans=imps[0].getNChannels();
				nslices=imps[0].getNSlices();
			}
			for(int j=0;j<stacks[0].getSize();j++){
				Object[] pix=new Object[ncols];
				for(int k=0;k<ncols;k++){
					pix[k]=jsmooth.bin2D(stacks[k].getPixels(j+1),width,height,binby,true);
				}
				Object temp=makeMontage(pix,width/binby,height/binby);
				stack.addSlice(name,temp);
			}
			IJ.showProgress(i,listtable.size());
		}
		jutils.create_hyperstack("Montage_Stack",stack,listtable.size(),nslices,nchans,true,null).show();
	}

	public Object makeMontage(Object[] imgs,int width,int height){
		if(imgs[0] instanceof byte[]){
			int newwidth=width*imgs.length;
			byte[] newimg=new byte[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((byte[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		} else if(imgs[0] instanceof short[]){
			int newwidth=width*imgs.length;
			short[] newimg=new short[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((short[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		} else {
			int newwidth=width*imgs.length;
			float[] newimg=new float[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((float[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		}
	}

}
