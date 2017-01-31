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
import ij.plugin.frame.*;
import ij.plugin.filter.*;

public class flatten_shifted_labels_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		RoiManager rman=RoiManager.getInstance();
		//rman.runCommand("show none");
		rman.runCommand("show all without labels");
		IJ.wait(100);
		ImagePlus imp2=imp.flatten();
		//(imp.flatten()).show();
		//IJ.selectWindow(imp.getTitle());
		IJ.wait(100);
		//rman.runCommand("show all");
		rman.runCommand("show all with labels");
		Roi[] rois=rman.getRoisAsArray();
		Filler filler=new Filler();
		for(int i=0;i<rois.length;i++){
			Rectangle r=new Rectangle(rois[i].getBounds().x-10,rois[i].getBounds().y-10,20,20);
			filler.drawLabel(imp2,imp2.getProcessor(),i+1,r);
		}
		imp2.show();
	}

}
