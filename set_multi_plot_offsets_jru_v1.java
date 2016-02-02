/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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

public class set_multi_plot_offsets_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("X Offset List");
		gd.addTextAreas("",null,10,20);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String input=gd.getNextText();
		String[] list=(new delimit_string(' ')).getrows(input);
		for(int i=0;i<npts.length;i++){
			float offset=Float.parseFloat(list[i]);
			for(int j=0;j<npts[i];j++){
				xvals[i][j]-=offset;
			}
		}
		jutils.runPW4VoidMethod(iw,"updatePlot");
	}

}
