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
import java.awt.Frame;
import java.util.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.io.*;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class stitch_pe_library_jru_v2 implements PlugIn,gui_interface {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		LOCI_file_reader r=new LOCI_file_reader();
		int nseries=r.getNSeries(directory,fname);
		String[] names=r.getSeriesNames(directory,fname);
		String[] keys={"X Location","Y Location"};
		String[][] vals=r.batch_get_series_metadata_value(directory,fname,keys);
		//now we group names into tiled sets (images end with '(raw tile x--)')
		String[] parentnames=new String[names.length];
		int[] tilenum=new int[names.length];
		float[][] icoords=new float[2][names.length];
		IJ.log("series name dump");
		for(int i=0;i<names.length;i++){
			int pos=names[i].indexOf("(raw tile ");
			if(pos<0){
				parentnames[i]="";
			}else{
				parentnames[i]=names[i].substring(0,pos-1);
				int pos2=names[i].indexOf(")",pos+8);
				String temp=names[i].substring(pos+10,pos2);
				tilenum[i]=Integer.parseInt(temp);
				icoords[0][i]=Float.parseFloat(vals[0][i]);
				icoords[1][i]=Float.parseFloat(vals[1][i]);
			}
			IJ.log(parentnames[i]+" , "+tilenum[i]+" , "+icoords[0][i]+" , "+icoords[1][i]);
		}
		//now get the unique names
		Object[] temptilenames=getUniqueNames(parentnames);
		String[] tilenames=(String[])temptilenames[0];
		String[] disptilenames=(String[])temptilenames[2];
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Select Image",disptilenames,disptilenames[0]);
		//gd.addNumericField("Overlap Percent",10.0,5,15,null);
		//gd.addNumericField("X Images",12,0);
		//gd.addNumericField("Y Images",8,0);
		gd.addCheckbox("Coordinates_from_plot",false);
		gd.addNumericField("Bin By",1,0);
		gd.addCheckbox("Subtract Avg Back",true);
		gd.addCheckbox("Smooth Back",false);
		gd.addNumericField("Smooth_Stdev",100,5,15,null);
		gd.addCheckbox("Output_Unstitched",false);
		gd.addCheckbox("Output_Coords",false);
		//String[] startoptions={"Upper_Left","Upper_Right","Lower_Left","Lower_Right"};
		//gd.addChoice("Scan_Start",startoptions,startoptions[2]);
		gd.addCheckbox("Single_Channel",false);
		gd.addNumericField("Selected_Channnel",0,0);
		gd.addCheckbox("Z_Project",false);
		gd.addCheckbox("Open_all_frames",true);
		gd.addNumericField("Frames_to_open",1,0);
		gd.addNumericField("Start_Frame",0,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int selindex=gd.getNextChoiceIndex();
		//float overlap=(float)(gd.getNextNumber()/100.0);
		//int ximgs=(int)gd.getNextNumber();
		//int yimgs=(int)gd.getNextNumber();
		boolean plotwindow=gd.getNextBoolean();
		int binby=(int)gd.getNextNumber();
		boolean sub=gd.getNextBoolean();
		boolean smoothback=gd.getNextBoolean();
		float sbstdev=(float)gd.getNextNumber();
		boolean dontstitch=gd.getNextBoolean();
		boolean outcoords=gd.getNextBoolean();
		//int scandiroption=gd.getNextChoiceIndex();
		boolean singchan=gd.getNextBoolean();
		int selchan=(int)gd.getNextNumber();
		boolean zproj=gd.getNextBoolean();
		boolean allframes=gd.getNextBoolean();
		int importframes=(int)gd.getNextNumber();
		int startframe=(int)gd.getNextNumber();
		//int xstart=0; int ystart=yimgs-1; int xscandir=1; int yscandir=-1;
		//if(scandiroption==0){ystart=0; yscandir=1;}
		//if(scandiroption==1){xstart=ximgs-1; xscandir=-1; ystart=0; yscandir=1;}
		//if(scandiroption==3){xstart=ximgs-1; xscandir=-1;}
		//now that we have the parent name, make our series list
		int[] selseries=new int[names.length];
		int[] seltilenum=new int[names.length];
		int ntiles=0;
		for(int i=0;i<names.length;i++){
			if(parentnames[i].equals(tilenames[selindex])){
				selseries[ntiles]=i;
				seltilenum[ntiles]=tilenum[i];
				ntiles++;
			}
		}
		selseries=trunc_array(selseries,ntiles);
		seltilenum=trunc_array(seltilenum,ntiles);
		int[] order=jsort.get_javasort_order(seltilenum);
		
		float[][] coords=new float[2][ntiles];
		for(int i=0;i<ntiles;i++){
			coords[0][i]=icoords[0][selseries[order[i]]]; coords[1][i]=icoords[1][selseries[order[i]]];
			//coords[0][i]/=(float)binby;
			//coords[1][i]/=(float)binby;
		}

		//now read in the huge image stack
		Object[] hugestack=null;
		int width=0; int height=0; int nchan=0; float psize=0.0f; int stacksize=0; int nslices=0;
		int[] limits={0,-1,0,-1,0,-1};
		if(!allframes){limits[4]=startframe; limits[5]=startframe+importframes-1;}
		IJ.showStatus("reading and binning tiles");
		for(int i=0;i<ntiles;i++){
			ImagePlus imp=(new LOCI_file_reader()).get_loci_subimp(directory,fname,false,selseries[order[i]],false,"Max",-1,limits);
			if(singchan) imp=getchanimp(imp,selchan);
			if(binby>1) imp=binimp(imp,binby);
			if(zproj) imp=zprojimp(imp);
			ImageStack tstack=imp.getStack();
			if(hugestack==null){
				width=imp.getWidth(); height=imp.getHeight(); nchan=imp.getNChannels();
				psize=(float)jutils.get_psize(imp); stacksize=tstack.getSize();
				nslices=imp.getNSlices(); if(nslices==1) nslices=imp.getNFrames();
				hugestack=new Object[stacksize*ntiles];
				//coords=getTileCoords(ximgs,yimgs,width,height,overlap,xstart,ystart,xscandir,yscandir);
				//coords=stitching.getTileCoords(ximgs,yimgs,width,height,overlap,xstart,ystart,xscandir,yscandir);
				if(outcoords){
					new PlotWindow4("Stitching_Coords","x","y",coords[0],coords[1]).draw();
					//return;
				}
			}
			for(int j=0;j<stacksize;j++) hugestack[i*stacksize+j]=tstack.getPixels(j+1);
			System.gc();
		}
		if(plotwindow){
			ImageWindow[] iw=jutils.selectPlots(false,1,new String[]{"coord_plot"});
			if(iw==null) return;
			float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw[0],"getXValues");
			float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
			int sel=(Integer)jutils.runPW4VoidMethod(iw[0],"getSelected"); if(sel<0) sel=0;
			for(int i=0;i<xvals[0].length;i++){
				coords[0][i]=xvals[sel][i]*psize;
				coords[1][i]=yvals[sel][i]*psize;
			}
		}
		//find the offsets
		float[] overlap1=stitching.getAvgOverlap(coords,width,height);
		IJ.log("avg overlap = "+overlap1[0]);
		float overlap=overlap1[0]/((float)width*psize);
		//now (optionally) calculate the average image
		if(sub){
			IJ.showStatus("calculating background");
			float[][] avgstack=getTAvg(hugestack,nchan,nslices,ntiles);
			if(smoothback){
				for(int i=0;i<avgstack.length;i++){
					jsmooth.blur2D(avgstack[i],sbstdev,width,height);
				}
			}
			//new ImagePlus("avg_stack",jutils.array2stack(avgstack,width,height)).show();
			//now subtract from the original hugeimp, dismantling it to save space as we go
			IJ.showStatus("subtracting the background");
			int counter=0;
			for(int i=0;i<ntiles;i++){
				for(int j=0;j<nslices;j++){
					for(int k=0;k<nchan;k++){
						float[] subpix=avgstack[k+j*nchan];
						float[] orpix=algutils.convert_arr_float(hugestack[counter]);
						for(int l=0;l<subpix.length;l++) orpix[l]-=subpix[l];
						hugestack[counter]=orpix;
						System.gc();
						counter++;
					}
				}
			}
			//new ImagePlus("sub_stack",jutils.array2stack(hugestack,width,height)).show();
		}
		String outname=tilenames[selindex].replace('\\','_');
		outname=outname.replace('/','_');
		if(dontstitch){
			ImagePlus rawimp=new ImagePlus(outname+"_raw",jutils.array2stack(hugestack,width,height));
			rawimp.setOpenAsHyperStack(true);
			rawimp.setDimensions(nchan,nslices,ntiles);
			jutils.set_psize(rawimp,psize);
			rawimp.show();
		} else {
			//now stitch it
			int typeindex=algutils.get_array_type(hugestack[0]);
			IJ.showStatus("stitching the image");
			stitching sclass=new stitching(width,height,coords[0],coords[1],psize,typeindex,this);
			ImageStack stitchstack=new ImageStack(sclass.newwidth,sclass.newheight);
			for(int i=0;i<nslices;i++){
				for(int j=0;j<nchan;j++){
					Object[] tser=algutils.get3DTSeries(hugestack,i,j,ntiles,nslices,nchan);
					stitchstack.addSlice("",sclass.stitch_frame(tser,true,overlap));
				}
			}
			hugestack=null;
			System.gc();

			ImagePlus stitchimp=new ImagePlus(outname+"_stitched",stitchstack);
			stitchimp.setOpenAsHyperStack(true);
			stitchimp.setDimensions(nchan,nslices,1);
			jutils.set_psize(stitchimp,psize);
			stitchimp.show();
		}
		IJ.showStatus("done stitching");
	}

	public static String getStitchedNames(String directory,String fname){
		int nseries=(new LOCI_file_reader()).getNSeries(directory,fname);
		String[] names=(new LOCI_file_reader()).getSeriesNames(directory,fname);
		//now we group names into tiled sets (images end with '(raw tile x--)')
		String[] parentnames=new String[names.length];
		int[] tilenum=new int[names.length];
		for(int i=0;i<names.length;i++){
			int pos=names[i].indexOf("(raw tile ");
			if(pos<0){
				parentnames[i]="";
			}else{
				parentnames[i]=names[i].substring(0,pos-1);
				int pos2=names[i].indexOf(")",pos+8);
				String temp=names[i].substring(pos+10,pos2);
				tilenum[i]=Integer.parseInt(temp);
			}
			//IJ.log(parentnames[i]+" , "+tilenum[i]);
		}
		//now get the unique names
		Object[] tilenames=getUniqueNames(parentnames);
		return table_tools.print_string_array((String[])tilenames[2]);
		//return tilenames;
	}

	public ImagePlus binimp(ImagePlus imp,int binby){
		ImageStack stack=imp.getStack();
		float psize=(float)jutils.get_psize(imp);
		int width=stack.getWidth(); int height=stack.getHeight();
		int slices=stack.getSize();
		ImageStack retstack=new ImageStack(width/binby,height/binby);
		for(int i=0;i<slices;i++){
			retstack.addSlice("",jsmooth.bin2D(stack.getPixels(i+1),width,height,binby,true));
		}
		ImagePlus retimp=new ImagePlus("binned",retstack);
		retimp.setOpenAsHyperStack(true);
		retimp.setDimensions(imp.getNChannels(),imp.getNSlices(),imp.getNFrames());
		retimp.copyScale(imp);
		jutils.set_psize(retimp,psize*(float)binby);
		return retimp;
	}

	public ImagePlus zprojimp(ImagePlus imp){
		ImageStack stack=imp.getStack();
		int width=stack.getWidth(); int height=stack.getHeight();
		int nchan=imp.getNChannels(); int nslices=imp.getNSlices(); int nframes=imp.getNFrames();
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nchan;j++){
				float[] proj=jutils.get3DProjZStat(stack,i,j,nframes,nslices,nchan,"Max");
				retstack.addSlice("",proj);
			}
		}
		ImagePlus retimp=new ImagePlus("Z Proj",retstack);
		retimp.setOpenAsHyperStack(true);
		retimp.setDimensions(nchan,1,nframes);
		retimp.copyScale(imp);
		return retimp;
	}

	public ImagePlus getchanimp(ImagePlus imp,int selchan){
		ImageStack stack=imp.getStack();
		int nchan=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		int counter=selchan;
		ImageStack retstack=new ImageStack(stack.getWidth(),stack.getHeight());
		for(int i=0;i<nframes*nslices;i++){
			retstack.addSlice("",stack.getPixels(counter+1));
			counter+=nchan;
		}
		ImagePlus retimp=new ImagePlus("singchan",retstack);
		retimp.setOpenAsHyperStack(true);
		retimp.setDimensions(1,nslices,nframes);
		retimp.copyScale(imp);
		return retimp;
	}

	public float[][] getTAvg(Object[] stack,int nchan,int nslices,int nframes){
		float[][] avg=new float[nchan*nslices][];
		int counter=0;
		for(int i=0;i<nslices;i++){
			for(int j=0;j<nchan;j++){
				avg[counter]=algutils.get3DProjTStat(stack,i,j,nframes,nslices,nchan,"Avg",null);
				counter++;
			}
		}
		return avg;
	}

	public int[] trunc_array(int[] arr,int len){
		int[] temp=new int[len];
		System.arraycopy(arr,0,temp,0,len);
		return temp;
	}

	public static Object[] getUniqueNames(String[] list1){
		String[] list=list1.clone();
		//sort the list
		jsort.javasort_order(list);
		List<String> celllist=new ArrayList<String>();
		//first skip all of the file names that aren't tiles (set to null)
		int pos=0;
		while(list[pos].length()==0){
			pos++;
		}
		String currcell=list[pos];
		celllist.add(currcell);
		int[] counts=new int[list1.length];
		counts[0]=1;
		for(int i=pos+1;i<list.length;i++){
			if(!list[i].equals(currcell)){
				celllist.add(list[i]);
				currcell=list[i];
			}
			counts[celllist.size()-1]++;
		}
		String[] temp=new String[celllist.size()];
		int[] tempcounts=new int[celllist.size()];
		String[] temp2=new String[celllist.size()];
		for(int i=0;i<celllist.size();i++){
			temp[i]=celllist.get(i);
			tempcounts[i]=counts[i];
			temp2[i]=temp[i]+" ("+tempcounts[i]+" tiles)";
		}
		return new Object[]{temp,tempcounts,temp2};
	}

	public void showMessage(String message){ IJ.showMessage(message);}

	public void showProgress(int currpos,int finalpos){ IJ.showProgress(currpos,finalpos);}

}
