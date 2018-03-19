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
import ij.io.*;
import java.io.*;
import jalgs.*;
import jguis.*;
import ij.text.*;

public class import_flowcyte_jru_v1 implements PlugIn {

	public void run(String arg) {
		//start by reading the header
		jdataio jdio=new jdataio();
		OpenDialog od = new OpenDialog("Open File",arg);
		String directory=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Show Metadata?",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean showmetadata=gd.getNextBoolean();
		try{
			//first read the header
			InputStream instream=new BufferedInputStream(new FileInputStream(directory+name));
			String label=jdio.readstring(instream,6);
			jdio.skipstreambytes(instream,4);
			String toff=jdio.readstring(instream,8); int textoff=Integer.parseInt(toff.trim()); //start of TEXT segment
			String teoff=jdio.readstring(instream,8); int texteoff=Integer.parseInt(teoff.trim()); //end of TEXT segment
			String doff=jdio.readstring(instream,8); int dataoff=Integer.parseInt(doff.trim()); //start of DATA segment
			String deoff=jdio.readstring(instream,8); int dataeoff=Integer.parseInt(deoff.trim()); //end of DATA segment
			String aoff=jdio.readstring(instream,8); int analoff=Integer.parseInt(aoff.trim()); //start of ANALYSIS segment
			String aeoff=jdio.readstring(instream,8); int analeoff=Integer.parseInt(aeoff.trim()); //end of ANALYSIS segment
			//IJ.log(""+textoff+" , "+texteoff+" , "+dataoff+" , "+dataeoff+" , "+analoff+" , "+analeoff);
			instream.close();
			instream=new BufferedInputStream(new FileInputStream(directory+name));
			jdio.skipstreambytes(instream,textoff);
			//now read all parameters and values
			int textlength=texteoff-textoff;
			String text=jdio.readstring(instream,textlength);
			jdio.skipstreambytes(instream,1);
			text=text.replace('\\','/');
			text=text.replace('|','/');
			text=text.replace("\u000C","/");
			text=text.replace("\u001E","/");
			String[] params=text.split("/\\u0024");
			//String[] params=text.split("/\u0c24");
			//IJ.log(text);
			String[] labels=new String[params.length-1];
			String[] values=new String[params.length-1];
			for(int i=1;i<params.length;i++){
				String[] temp=params[i].split("/");
				labels[i-1]=temp[0]; values[i-1]=temp[1];
				if(showmetadata){
					IJ.log(labels[i-1]+" , "+values[i-1]);
				}
			}
			//now get the number of data points and number of channels
			int npts=get_label_number(labels,values,"TOT");
			int nch=get_label_number(labels,values,"PAR");
			String byteord=get_label_value(labels,values,"BYTEORD");
			//IJ.log(byteord);
			boolean motorola=true;
			if(byteord.startsWith("1,")) motorola=false;
			//now go through each channel getting its name, bits, and range
			String[] chnames=new String[nch];
			int[] chbits=new int[nch];
			int[] range=new int[nch];
			for(int i=0;i<nch;i++){
				chnames[i]=get_label_value(labels,values,"P"+(i+1)+"N");
				chbits[i]=get_label_number(labels,values,"P"+(i+1)+"B");
				range[i]=get_label_number(labels,values,"P"+(i+1)+"R");
			}
			String headings=table_tools.print_string_array(chnames);
			TextWindow tw=new TextWindow(name,headings,"",400,200);
			float[][] temp=new float[npts][nch];
			instream.close();
			instream=new BufferedInputStream(new FileInputStream(directory+name));
			jdio.skipstreambytes(instream,dataoff);
			for(int i=0;i<npts;i++){
				for(int j=0;j<nch;j++){
					if(chbits[j]==16){
						if(!motorola) temp[i][j]=(float)jdio.readintelshort(instream);
						else temp[i][j]=(float)jdio.readmotorolashort(instream);
					} else {
						if(chbits[j]==32){
							if(!motorola) temp[i][j]=(float)jdio.readintelfloat(instream);
							else temp[i][j]=(float)jdio.readmotorolafloat(instream);
						} else {
							if(chbits[j]==8){
								temp[i][j]=(float)jdio.readintelbyte(instream);
							}
						}
					}
				}
			}
			String data=table_tools.print_float_array(temp);
			tw.append(data);
			instream.close();
		} catch(IOException e){
			IJ.log(e.getMessage());
		}
	}

	private String get_label_value(String[] labels,String[] values,String match){
		for(int i=0;i<labels.length;i++){
			if(labels[i].equals(match)){
				return values[i];
			}
		}
		return null;
	}

	private int get_label_number(String[] labels,String[] values,String match){
		String value=get_label_value(labels,values,match);
		if(value!=null){
			try{
				return Integer.parseInt(value.trim());
			}catch(NumberFormatException e){
				return -1;
			}
		}
		return -1;
	}

}
