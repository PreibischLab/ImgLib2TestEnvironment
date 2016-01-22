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
	public static <T extends RealType<T>> PointSampleList<T> getsubList(PointSampleList<T> list, int direction) {

		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
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

		final Cursor<T> listCursor = list.localizingCursor().copyCursor();

		PointSampleList<T> childA = new PointSampleList<T>(n);

		PointSampleList<T> childB = new PointSampleList<T>(n);

		PointSampleList<T> sortedList = new PointSampleList<T>(n);

		meanIndex = (int) (list.min(direction) + ((list.max(direction) - list.min(direction)) / 2));

		// In this loop I create the splitPlane at mean value in one direction

		while (listCursor.hasNext()) {

			listCursor.fwd();

			Point splitPoint = new Point(n);

			long splitPlane = listCursor.getLongPosition(otherdirection);

			splitPoint.setPosition(meanIndex, direction);
			splitPoint.setPosition(splitPlane, otherdirection);

			Point cord = new Point(n);

			cord.setPosition(listCursor);
			if (listCursor.getLongPosition(direction) <= splitPoint.getDoublePosition(direction))

				childA.add(cord, listCursor.get().copy());

			else

				childB.add(cord, listCursor.get().copy());

		}

		if ((meanIndex) >= 0 && childA.size() > 0 && childB.size() > 0) {

			// System.out.println(" Size of list childA: "+childA.size());
			// System.out.println(" Size of list childB: "+childB.size());

			getsubList(childA, otherdirection);
			getsubList(childB, otherdirection);

			sortedList = mergeList(childA, childB);
		}

		return sortedList;

	}

	public static <T extends RealType<T>> PointSampleList<T> mergeList(PointSampleList<T> listA,
			PointSampleList<T> listB) {

		int n = listA.numDimensions();
		int m = listB.numDimensions();

		PointSampleList<T> mergedList = new PointSampleList<T>(n);

		Cursor<T> cursorA = listA.cursor();

		Cursor<T> cursorB = listB.cursor();

		cursorA.fwd();
		cursorB.fwd();

		if (n == m) {

			while (cursorA.hasNext()) {

				Point cord = new Point(n);

				cord.setPosition(cursorA);

				cursorA.fwd();
				cursorB.fwd();
				if (cursorA.get().compareTo(cursorB.get()) > 0) {

					cursorA.get().set(cursorB.get().copy());

					mergedList.add(cord, cursorA.get().copy());
					mergedList.add(cord, cursorB.get().copy());
					
				} else {

					cursorA.get().copy();

					mergedList.add(cord, cursorA.get().copy());
					mergedList.add(cord, cursorB.get().copy());
					

				}
				
				
			}
		}
		
		
		
	

		return mergedList;

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

		PointSampleList<FloatType> sortedlist = new PointSampleList<FloatType>(img.numDimensions());
		
		
		

		list = getList(img);

		sortedlist = getsubList(list, 0);

		
		
		
	Cursor<FloatType> test= sortedlist.cursor();
	Cursor<FloatType> testiter= sortedlist.cursor();
		
		test.fwd();
		testiter.fwd();
		
		while(test.hasNext()){
		testiter.next();
		
		if(test.get().compareTo(testiter.get())>0 )
			
			System.out.println("False");
		
		
		
		test=testiter;
		}
			
			
			
				
				
		
		
		System.out.println(list.size());
		System.out.println(sortedlist.size());

		// ImageJFunctions.show(img).setTitle("Original_Image");

		// ImageJFunctions.show(imgout).setTitle("Permuted_Image");

	}

}
