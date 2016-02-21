package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.neighborhood.DiamondTipsNeighborhood.LocalCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;

import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class distanceTransform {

	public static void createBitimage(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<BitType> imgout,
			FloatType ThresholdValue) {

		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<BitType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);
		
			if (bound.get().compareTo(ThresholdValue) > 0) {

				outbound.get().setOne();

			}

			else {

				outbound.get().setZero();

			}

		}
	}

	
	public interface Distance{
		
		
		double getdistance(Localizable cursor1, Localizable cursor2);
		
		
	}
	
	public static class EucledianDistance implements Distance
	{
	public  double getdistance(Localizable cursor1, Localizable cursor2){
		
		double distance=0;
		
		for (int d=0; d<cursor2.numDimensions(); ++d){
			
			distance += Math.pow((cursor2.getDoublePosition(d) - cursor1.getDoublePosition(d)), 2);
		}
		
		
		
		return Math.sqrt(distance);
		
		
		
		
	}
	
	}
	
	
	public static class MannhattanDistance implements Distance
	{
	public  double getdistance(Localizable cursor1, Localizable cursor2){
		
		double distance=0;
		
		for (int d=0; d<cursor2.numDimensions(); ++d){
			
			distance += Math.abs(cursor2.getDoublePosition(d) - cursor1.getDoublePosition(d));
		}
		
		
		
		return distance;
		
		
		
		
	}
	
	}
	
	
	public static <T extends RealType<T>> void computeDistance(RandomAccessibleInterval<BitType> img,
			RandomAccessibleInterval<T> imgout, final Distance dist) throws FileNotFoundException {

		final Cursor<BitType> bound = Views.iterable(img).cursor();

		final RandomAccess<T> outbound = imgout.randomAccess();
		PrintStream out = new PrintStream(new FileOutputStream("BruteForcedist.txt"));
		System.setOut(out);
		while (bound.hasNext()) {
			bound.fwd();
			outbound.setPosition(bound);

		
			if (bound.get().getInteger() == 0) {
				final Cursor<BitType> second = Views.iterable(img).cursor();
				double mindistance = Double.MAX_VALUE;
				

				double distance = 0;
				while (second.hasNext()) {
					if (second.next().getInteger() == 1) {

						distance=dist.getdistance(bound,second);

				mindistance = Math.min(mindistance, distance);
				
					}
				}
				System.out.println(mindistance);

				outbound.get().setReal(mindistance);
			
				
				}

			else {
				outbound.get().setReal(0);
			}

		}

	}

	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));

		ImageJFunctions.show(img).setTitle("Original_Image");

		final Img<BitType> imgout = new ArrayImgFactory<BitType>().create(img, new BitType());

		final Img<BitType> bitimgout = new ArrayImgFactory<BitType>().create(img, new BitType());

		FloatType val = new FloatType(200);

		createBitimage(img, imgout, val);

		computeDistance(imgout, img, new EucledianDistance());

		ImageJFunctions.show(img).setTitle("Eucledian_FloatType_output");
/*
		computeDistance(imgout, bitimgout, new EucledianDistance());

		ImageJFunctions.show(bitimgout).setTitle("Eucledian_BitType_output");
		
		
		computeDistance(imgout, img, new MannhattanDistance());

		ImageJFunctions.show(img).setTitle("Mannhattan_FloatType_output");
		
		computeDistance(imgout, bitimgout, new MannhattanDistance());

		ImageJFunctions.show(bitimgout).setTitle("Mannhattan_BitType_output");
		*/
		

	}

}
