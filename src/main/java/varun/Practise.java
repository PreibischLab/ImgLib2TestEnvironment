package varun;

import java.util.Random;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.Sampler;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class Practise {
	
	public static < T extends RealType < T > > void newSphere
	(RandomAccessibleInterval < T > img, Point center, long radius) {
	
		
		HyperSphere<T> mysphere = new HyperSphere<T> (img, center, radius);
		
		long defaultcenter[]=new long[img.numDimensions()];
		
		for(int d=0;d<img.numDimensions();++d){
			
			defaultcenter[d]=(img.max(d)-img.min(d))/2+img.min(d);
			
		}
		Cursor<T> bigc = mysphere.cursor();
		while(bigc.hasNext()){
			bigc.fwd();
		}
		Random ranval= new Random();
		double ranlocation=ranval.nextDouble();
		long smallradius=Math.round(ranval.nextDouble());
		
		HyperSphere<T> smallSphere=new HyperSphere<T>(img, bigc, smallradius);
		
		Cursor < T > smallc =smallSphere.cursor();
		
		
		
		while(smallc.hasNext()){
			smallc.fwd();
			
//T value =   smallc.get();
//value.setReal(ranlocation);
smallc.get();
    
			
		}
		
	}

	public static void main (String[] args){
		
		final Img<FloatType> img = ArrayImgs.floats(1000, 1000);
	long center[]=new long[img.numDimensions()];
		
		for(int d=0;d<img.numDimensions();++d){
			
			center[d]=(img.max(d)-img.min(d))/2+img.min(d);
			
		}
		Point centerdef = new Point(center);
		
		long radius=300;
		newSphere(img,centerdef,radius);
		
	}
	
	
}
