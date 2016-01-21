package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import net.imglib2.PointSampleList;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.kdtree.SplitHyperPlaneKDTree;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Kdbranch {

	

	// Input a List and get two sublists split at the median
	public static <T extends RealType<T>> void getsubList(PointSampleList<T> list, int direction) {

		
		/****  To ward against running over the dimensionality, creating some local restrictions on the global variable direction   ****/
		if (direction == list.numDimensions())
			direction = 0;
		
		int otherdirection;
		
		otherdirection = direction + 1;
		if (direction == list.numDimensions() - 1)
			otherdirection = 0;
		if (direction >= 0 && direction < list.numDimensions() - 1)
			otherdirection = direction + 1;
		
		int n = list.numDimensions();

		int meanIndex;
		

		final Cursor<T> listCursorA = list.localizingCursor().copyCursor();
		final Cursor<T> listCursorB = list.localizingCursor().copyCursor();
		
		

		PointSampleList<T> childA = new PointSampleList<T>(n);

		PointSampleList<T> childB = new PointSampleList<T>(n);

		meanIndex = (int) ((list.max(direction) - list.min(direction)) / 2);
		
		//System.out.println("meanIndex: "+meanIndex);
		
		// In this loop I create the splitPlane at mean value in one direction, this is for childA which always has the point 0,0...
		
		while (listCursorA.hasNext()) {

			listCursorA.fwd();
			
			Point splitPoint = new Point(n);

			long splitPlane = listCursorA.getLongPosition(otherdirection);

			splitPoint.setPosition(meanIndex, direction);
			splitPoint.setPosition(splitPlane, otherdirection);
			
			Point cord = new Point(n);

			cord.setPosition(listCursorA);
			if (listCursorA.getLongPosition(direction) <= splitPoint.getDoublePosition(direction))
				childA.add(cord, listCursorA.get().copy());
			
			
		}
		
		
		// For sake of clarity I create a separate loop to make the childB which does-not have the point 0,0...
		
		
		while (listCursorB.hasNext()) {

			listCursorB.fwd();
			
			Point splitPoint = new Point(n);

			long splitPlane = listCursorB.getLongPosition(otherdirection);

			splitPoint.setPosition(meanIndex, direction);
			splitPoint.setPosition(splitPlane, otherdirection);
			
			Point cord = new Point(n);

			cord.setPosition(listCursorB);
			
			if (listCursorB.getLongPosition(direction) > splitPoint.getDoublePosition(direction))
				childB.add(cord, listCursorB.get().copy());
			
			
		}
		
		//This is the partition having the point 0,0.. This recursion works fine
		if ((meanIndex - 1) > 0 && childA.size() > 0) {
			System.out.print("Number of pixels in current direction: "+list.dimension(direction));
			
			System.out.println("    meanIndex: "+meanIndex);
			
			getsubList(childA, otherdirection);
			
		}
		
		
	/***** Problem part, Can not create right side of the tree which does not have the point 0,0.. ******/  
		
		// This is the partition not having the point 0,0... and I think for this reason the space is not further partitioned.
		
		
		/*
		if ((meanIndex + 1) <= (list.dimension(direction) - 1) && childB.size() > 0) {
			 getsubList(childB, otherdirection); 
			 
		}
*/
		
	}

	public static <T extends RealType<T>> PointSampleList<T> getList(RandomAccessibleInterval<T> img) {

		final RandomAccessible<T> infinite = Views.extendZero(img);

		final int n = img.numDimensions();
		long min[] = new long[n];
		long max[] = new long[n];

		for (int d = 0; d < n; ++d) {

			min[d] = img.min(d);
			max[d] = img.max(d);

		}

		FinalInterval interval = new FinalInterval(min, max);

		final IterableInterval<T> imgav = Views.interval(infinite, interval);

		final Cursor<T> first = imgav.cursor();

		// A point sample list with coordinates declared and initialized.
		PointSampleList<T> parent = new PointSampleList<T>(n);

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(n);

			cord.setPosition(first);

			parent.add(cord, first.get().copy());

		}

		return parent;

	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		list = getList(img);

		getsubList(list, 0);

		// ImageJFunctions.show(img).setTitle("Original_Image");

		// ImageJFunctions.show(imgout).setTitle("Permuted_Image");

	}

}
