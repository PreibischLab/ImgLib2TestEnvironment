package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class Ex6a1 {
	
	public Ex6a1(){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		double[] sigma = new double[img.numDimensions()];
		for (int i = 2; i <= 8; i*=2){
			for (int d = 0; d < img.numDimensions(); ++d)
				sigma[d] = i;
			ImageJFunctions.show(Gauss.toFloat(sigma, img));
		}
	}
	
	public static void main(String[] args){
		new ImageJ();
		new Ex6a1();
	}
}
