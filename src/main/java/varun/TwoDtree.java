package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Random;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.ui.util.StopWatch;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class TwoDtree {

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

	// Sorts the co-ordinates in a given direction, the central element is then
	// always the pivot for the kDTree.
	public static <T extends RealType<T>> ArrayList<Long> sortedCoordinates(PointSampleList<T> list, int direction) {

		final Cursor<T> listCursor = list.localizingCursor();

		final ArrayList<Long> values = new ArrayList<Long>((int) list.dimension(direction));

		while (listCursor.hasNext()) {
			listCursor.next();
			values.add(listCursor.getLongPosition(direction));
		}

		split(values, direction);

		return values;

	}

	public static <T extends RealType<T>> void splitbyCoordinate(PointSampleList<T> list,
			ArrayList<Long> sortedcoordinateList,int startindex, int lastindex, int direction) {

		int n = list.numDimensions();
		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == list.numDimensions())
			direction = 0;
		if (list.dimension(direction) <= 1)
			return;

		else {

			/****
			 * To ward against running over the dimensionality, creating some
			 * local restrictions on the global variable direction
			 ****/
			if (direction == list.numDimensions())
				direction = 0;

			// the first element belonging to the right list childB

			final PointSampleList<T> childA = new PointSampleList<T>(n);
			final PointSampleList<T> childB = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			double pivot;
			
			



			int startindexA= 0;
			int lastindexA=(int) (list.dimension(direction)/2-1+list.dimension(direction)%2);
			
			int startindexB=(int) (list.dimension(direction)/2+list.dimension(direction)%2);
			int lastindexB=(int) list.dimension(direction)-1;

			

			 pivot= pivotElement(sortedcoordinateList,
			 direction,startindex,lastindex);

			while (listCursor.hasNext()) {

				listCursor.fwd();
				

				Point cord = new Point(listCursor);

				if (listCursor.getDoublePosition(direction) < pivot) {

					childA.add(cord, listCursor.get().copy());

					// System.out.println("childA: "+listCursor.get());

				} else

				{

					childB.add(cord, listCursor.get().copy());
					// System.out.println("childB: "+listCursor.get());
				}

			}

			// splitbyCoordinate(childA,startindexA,lastindexA, direction+1);

			// splitbyCoordinate(childB,startindexB, lastindexB, direction+1);

		}

	}
	
	
	public static <T extends RealType<T>> Localizable firstLocation(final IterableInterval<T> interval) {
		Cursor<T> c = interval.localizingCursor();
		c.fwd();
		
		return c;
	}
	

	public static <T extends RealType<T>> double pivotElement(ArrayList<Long> sortedcoordinateList, int direction,
			int startindex, int lastindex) {

		double Element;

		int pivotIndex = startindex + (lastindex - startindex) / 2;

		Element = sortedcoordinateList.get(pivotIndex);

		// System.out.println(Element);

		return Element;

	}

	public static <T extends RealType<T>> void split(ArrayList<Long> coordinateList, int direction) {

		if (coordinateList.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) coordinateList.size() / 2;

			Iterator<Long> iterator = coordinateList.iterator();

			final ArrayList<Long> childA = new ArrayList<Long>((int) coordinateList.size() / 2);

			final ArrayList<Long> childB = new ArrayList<Long>(
					(int) coordinateList.size() / 2 + coordinateList.size() % 2);

			int xindex = 0;

			while (iterator.hasNext()) {
				iterator.next();
				if (xindex < splitIndex) {
					childA.add(xindex, coordinateList.get(xindex));

				} else

					childB.add(xindex, coordinateList.get(xindex));

				xindex++;

			}

			// System.out.println("childA : " + childA.size());

			// System.out.println("childB : " + childB.size());

			split(childA, direction);

			split(childB, direction);

			mergeListValue(coordinateList, childA, childB);

		}
		System.out.println("Sorted List : " + coordinateList);
	}

	/// ***** Returns a sorted list *********////
	public static <T extends RealType<T>> void mergeListValue(ArrayList<Long> sortedlist, ArrayList<Long> listA,
			ArrayList<Long> listB) {

		int i = 0, j = 0, k = 0;

		while (i < listA.size() && j < listB.size()) {

			if (listA.get(i) < listB.get(j)) {

				sortedlist.set(k, listA.get(i));

				++i;
				++k;
			}

			else {

				sortedlist.set(k, listB.get(j));

				++j;
				++k;

			}

		}

		while (i < listA.size()) {
			sortedlist.set(k, listA.get(i));
			++i;
			++k;

		}

		while (j < listB.size()) {
			sortedlist.set(k, listB.get(j));
			++j;
			++k;

		}

	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		ArrayList<Long> XcoordinatesSort = new ArrayList<Long>((int) list.dimension(0));
		ArrayList<Long> YcoordinatesSort = new ArrayList<Long>((int) list.dimension(1));

		// Make a 1D list along the X direction by setting an appropriate
		// interval on the image.

		IterableInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 5, 1 });

		final Cursor<FloatType> first = view.cursor();

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(img.numDimensions());

			cord.setPosition(first);

			list.add(cord, first.get().copy());
			// System.out.println("Set of x co-ordinates Initial List : " +
			// cord.getDoublePosition(0));
			// System.out.println("Set of y co-ordinates Initial List : " +
			// cord.getDoublePosition(1));
			// System.out.println("Values Initial list : " + first.get());

		}

		XcoordinatesSort = sortedCoordinates(list, 0);
		YcoordinatesSort = sortedCoordinates(list, 1);
		splitbyCoordinate(list,XcoordinatesSort,0,(int)list.dimension(0), 0); // Split list along X direction
		splitbyCoordinate(list,YcoordinatesSort,0,(int)list.dimension(1), 1); // Split list along Y direction

		Cursor<FloatType> testtwo = list.cursor();

		while (testtwo.hasNext()) {
			testtwo.fwd();
			Point newpoint = new Point(img.numDimensions());

			newpoint.setPosition(testtwo);

			// System.out.println("Set of x co-ordinates sorted List : " +
			// newpoint.getDoublePosition(0));
			// System.out.println("Set of y co-ordinates sorted List : " +
			// newpoint.getDoublePosition(1));
			// System.out.println("Values sorted list : " + testtwo.get());

		}

	}

}
