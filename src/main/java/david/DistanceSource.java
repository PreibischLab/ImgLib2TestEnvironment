package david;

import ij.ImageJ;

import java.io.File;

import util.ImgLib2Util;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

/**
 * a class representing a local minimum in one row/column of the distance map-in-progress
 * providing means of updating the distance map with distances growing from the minimum (in one dimension)
 * 
 * @author David
 *
 */

public class DistanceSource implements Comparable<DistanceSource> {
	
	public enum DistanceType{
		EUCLIDEAN,
		MANHATTAN
	}
	
	long pos;
	long dist;
	double initDist;
	double nextDist;
	boolean validUp;
	boolean validDown;
	
	public DistanceSource(long pos, double initDist) {
		this.pos = pos;
		this.dist = 1;
		this.initDist = initDist;
		this.nextDist = Math.sqrt(initDist*initDist + 1);
		this.validDown = true;
		this.validDown = true;
	}
	
	public boolean isValid() {
		return validUp || validDown;
	}
	
	public <T extends RealType<T>> void applyToImageAndGrow(RandomAccessibleInterval<T> line, DistanceType distanceType) {
		RandomAccess<T> ra = line.randomAccess();
		
		// update DOWN (towards lower index)
		// set validity to false if end of line is reached or no update to lower distance possible any more 
		if (pos - dist >= 0){
			ra.setPosition(pos-dist, 0);
			if (ra.get().getRealDouble() > nextDist){
				ra.get().setReal(nextDist);
			} else {
				validDown = false;
			}
		} else {
			validDown = false;
		}
		
		// update UP
		if (pos + dist < line.dimension(0)){
			ra.setPosition(pos+dist, 0);
			if (ra.get().getRealDouble() > nextDist){
				ra.get().setReal(nextDist);
			} else {
				validUp = false;
			}
		} else {
			validUp = false;
		}
		
		// update next distance
		dist++;
		if (distanceType == DistanceType.EUCLIDEAN) {
			nextDist = Math.sqrt(initDist*initDist + dist*dist);
		} else if (distanceType == DistanceType.MANHATTAN){
			nextDist = initDist + dist;
		}		
		
	}

	@Override
	public int compareTo(DistanceSource arg0) {
		return Double.compare(this.nextDist, arg0.nextDist);
	}
	
	public static <T extends RealType<T>> RandomAccessibleInterval<T> getLine(RandomAccessibleInterval<T> img, int d){
		
		RandomAccessibleInterval<T> res = img;
		
		for (int i = 0; i < img.numDimensions(); i++){
			if ( i != d){
				res = Views.hyperSlice(res, i, 0);
			}
		}
		
		return res;
		
	}
	
	public static <T extends RealType<T>> void calcDT(RandomAccessibleInterval<T> dt, DistanceType dstType){
		
		
		for (int d = 0; d < dt.numDimensions(); d ++){			
			for (int i = 0; i< dt.numDimensions(); i++){
				
				RandomAccessibleInterval<T> tImg = dt;
				if ( i != d){
					for (int j = 0; j < dt.dimension(i); j++){
						tImg = Views.hyperSlice(tImg, i, j);
					}
				}
				
			}
			
			
		}
		
		
	}
	
	public static void main(String[] args) {
		new ImageJ();
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(
				"src/main/resources/bridge.png"));
		ImageJFunctions.show(img);
		Img<BitType> thresholded = new ArrayImgFactory<BitType>().create(img,
				new BitType());
		Thresholding.threshold(img, thresholded, new FloatType(200));
		ImageJFunctions.show(thresholded);
		
		RandomAccessibleInterval<FloatType> a = Views.dropSingletonDimensions(getLine(img, 1));
		
		ImageJFunctions.show(Views.dropSingletonDimensions(getLine(img, 1)));
		
		
	}

}


