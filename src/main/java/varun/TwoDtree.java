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
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
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
			listCursor.fwd();

			values.add(listCursor.getLongPosition(direction));

		}
		// System.out.println(values);
		split(values, direction);
		// System.out.println(values);
		return values;

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

				if (xindex < splitIndex)
					childA.add(coordinateList.get(xindex));

				else

					childB.add(coordinateList.get(xindex));

				xindex++;

			}

			// System.out.println("childA : " + childA.size());

			// System.out.println("childB : " + childB.size());

			split(childA, direction);

			split(childB, direction);

			mergeListValue(coordinateList, childA, childB);

		}
		// System.out.println("Sorted List : " + coordinateList);
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

	public static <T extends RealType<T>> ArrayList<Double> medianValueLeft(ArrayList<Long> sortedcoordinateList,
			int startindex, int lastindex, int direction) {

		// Size of a list is lastindex-startindex+1.

		int medianIndexA, medianIndexB;
		if ((lastindex - startindex + 1) % 2 == 1) {
			medianIndexA = startindex + (lastindex - startindex + 1 + (lastindex - startindex + 1) % 2) / 2;
			medianIndexB = medianIndexA;
		}

		else {

			medianIndexA = startindex + (lastindex - startindex + 1 + (lastindex - startindex + 1) % 2) / 2 - 1;
			medianIndexB = medianIndexA + 1;
		}
		double medianElement;
		ArrayList<Double> leftmedians = new ArrayList<Double>();

		medianElement = 0.5 * (sortedcoordinateList.get(medianIndexA) + sortedcoordinateList.get(medianIndexB));
		leftmedians.add(medianElement);

		while (lastindex > startindex) {

			int startindexleft = startindex;
			int lastindexleft = medianIndexA - 1;
			int medianIndexleftA, medianIndexleftB;

			double firstElement;
			firstElement = sortedcoordinateList.get(startindex);
			if (lastindexleft - startindexleft + 1 <= 2)
				break;

			if ((lastindexleft - startindexleft + 1) % 2 == 1) {
				medianIndexleftA = startindexleft
						+ (lastindexleft - startindexleft + 1 + (lastindexleft - startindexleft + 1) % 2) / 2;
				medianIndexleftB = medianIndexleftA;
			}

			else {

				medianIndexleftA = startindexleft
						+ (lastindexleft - startindexleft + 1 + (lastindexleft - startindexleft + 1) % 2) / 2 - 1;
				medianIndexleftB = medianIndexleftA + 1;

			}

			medianElement = 0.5
					* (sortedcoordinateList.get(medianIndexleftA) + sortedcoordinateList.get(medianIndexleftB));

			lastindex = lastindexleft;
			medianIndexA = medianIndexleftA;
			if (medianElement > firstElement)
				leftmedians.add(medianElement);
			else
				break;

		}

		return leftmedians;

	}

	public static <T extends RealType<T>> ArrayList<Double> medianValueRight(ArrayList<Long> sortedcoordinateList,
			int startindex, int lastindex, int direction) {

		int medianIndexA, medianIndexB;
		if ((lastindex - startindex + 1) % 2 == 1) {
			medianIndexA = startindex + (lastindex - startindex + 1 + (lastindex - startindex + 1) % 2) / 2;
			medianIndexB = medianIndexA;
		}

		else {

			medianIndexA = startindex + (lastindex - startindex + 1 + (lastindex - startindex + 1) % 2) / 2 - 1;
			medianIndexB = medianIndexA + 1;
		}

		double medianElement;
		ArrayList<Double> rightmedians = new ArrayList<Double>();

		medianElement = 0.5 * (sortedcoordinateList.get(medianIndexA) + sortedcoordinateList.get(medianIndexB));
		rightmedians.add(medianElement);

		while (lastindex > startindex) {

			int startindexright = medianIndexB + 1;
			int lastindexright = lastindex;
			int medianIndexrightA, medianIndexrightB;

			double lastElement;
			lastElement = sortedcoordinateList.get(lastindex);
			if (lastindexright - startindexright + 1 <= 2)
				break;

			if ((lastindexright - startindexright + 1) % 2 == 1) {
				medianIndexrightA = startindexright
						+ (lastindexright - startindexright + 1 + (lastindexright - startindexright + 1) % 2) / 2;
				medianIndexrightB = medianIndexrightA;
			}

			else {
				medianIndexrightA = startindexright
						+ (lastindexright - startindexright + 1 + (lastindexright - startindexright + 1) % 2) / 2 - 1;
				medianIndexrightB = medianIndexrightA + 1;

			}

			medianElement = 0.5
					* (sortedcoordinateList.get(medianIndexrightA) + sortedcoordinateList.get(medianIndexrightB));

			startindex = startindexright;
			medianIndexB = medianIndexrightB;

			if (lastElement > medianElement)
				rightmedians.add(medianElement);

			else
				break;

		}
		return rightmedians;

	}

	public static <T extends RealType<T>> PointSampleList<T> getLeftTree(PointSampleList<T> list,
			ArrayList<Double> medianElements, int direction) {
		int n = list.numDimensions();
		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == list.numDimensions())
			direction = 0;
		if (list.dimension(direction) <= 2)
			return null;

		else {
			double pivotElement;

			pivotElement = medianElements.get(0);
			final PointSampleList<T> LeftTree = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			// In the list of medianValues the starting index stores the root
			// node in the direction, proceeded by medianValues on the left side
			// of the tree.

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement)

					LeftTree.add(cord, listCursor.get().copy());

			}

			return LeftTree;

		}

	}

	public static <T extends RealType<T>> PointSampleList<T> getRightTree(PointSampleList<T> list,
			ArrayList<Double> medianElements, int direction) {
		int n = list.numDimensions();
		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == list.numDimensions())
			direction = 0;
		if (list.dimension(direction) <= 2)
			return null;

		else {
			double pivotElement;

			pivotElement = medianElements.get(0);
			final PointSampleList<T> RightTree = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			// In the list of medianValues the starting index stores the root
			// node in the direction, proceeded by medianValues on the left side
			// of the tree.

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

				if (listCursor.getDoublePosition(direction) >= pivotElement)

					RightTree.add(cord, listCursor.get().copy());

			}

			return RightTree;

		}

	}

	public static <T extends RealType<T>> Pair<PointSampleList<T>, PointSampleList<T>> getsubTrees(
			PointSampleList<T> LeftorRightTree, ArrayList<Double> medianElements, int medianindex, int direction) {

		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == LeftorRightTree.numDimensions())
			direction = 0;
		if (LeftorRightTree.dimension(direction) <= 2)
			return null;

		else {
			int n = LeftorRightTree.numDimensions();
			double pivotElement;
			final PointSampleList<T> childA = new PointSampleList<T>(n);
			final PointSampleList<T> childB = new PointSampleList<T>(n);

			final Cursor<T> listCursor = LeftorRightTree.localizingCursor();

			pivotElement = medianElements.get(medianindex);

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point newpoint = new Point(n);
				newpoint.setPosition(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement) {

					childA.add(newpoint, listCursor.get().copy());

					// System.out.println("Hi");

				} else

				{

					childB.add(newpoint, listCursor.get().copy());
					// System.out.println("Hello");
				}

			}

			// System.out.println(pivotElement);

			Pair<PointSampleList<T>, PointSampleList<T>> pair = new ValuePair<PointSampleList<T>, PointSampleList<T>>(
					childA, childB);

			return pair;

		}

	}

	public static <T extends RealType<T>> double[] searchTree(ArrayList<Double> medianElementsLeft,
			ArrayList<Double> medianElementRight, double coordinate) {

		double[] pointLocation = new double[2];
		if (coordinate < medianElementsLeft.get(0)) {
			System.out.println("The point is on the Left Tree");
			// The point is on the Left Tree
			for (int medianindex = 1; medianindex < medianElementsLeft.size(); ++medianindex) {

				// Keep Searching till you find the value it is greater or equal
				// to

				if (coordinate <= medianElementsLeft.get(medianindex)
						&& coordinate >= medianElementsLeft.get(medianindex - 1)) {
					// Now you store the two index values

					pointLocation[0] = medianElementsLeft.get(medianindex);
					pointLocation[1] = medianElementsLeft.get(medianindex - 1);

				}

				else {
					
					pointLocation[0]=medianElementsLeft.get(medianElementsLeft.size()-1);
					pointLocation[1]=medianElementsLeft.get(medianElementsLeft.size()-1);
					
				}
				
				
			}

		}

		else if (coordinate >= medianElementRight.get(0)) {
			System.out.println("The point is on the Right Tree");
			// The point is on the Right Tree
			for (int medianindex = 1; medianindex < medianElementRight.size(); ++medianindex) {

				// Keep Searching till you find the value it is greater or equal
				// to

				if (coordinate <= medianElementRight.get(medianindex)
						&& coordinate >= medianElementRight.get(medianindex - 1)) {
					// Now you store the two index values

					pointLocation[0] = medianElementRight.get(medianindex);
					pointLocation[1] = medianElementRight.get(medianindex - 1);

				}
				
				else{
					
					pointLocation[0]=medianElementRight.get(medianElementRight.size()-1);
					pointLocation[1]=medianElementRight.get(medianElementRight.size()-1);
					
				}

			}

		}

		return pointLocation;
	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		ArrayList<Long> XcoordinatesSort = new ArrayList<Long>((int) list.dimension(0));
		ArrayList<Long> YcoordinatesSort = new ArrayList<Long>((int) list.dimension(1));

		ArrayList<Double> MedianLeftX = new ArrayList<Double>();

		ArrayList<Double> MedianRightX = new ArrayList<Double>();

		ArrayList<Double> MedianLeftY = new ArrayList<Double>();

		ArrayList<Double> MedianRightY = new ArrayList<Double>();

		// Make a list by setting an appropriate
		// interval on the image.

		IterableInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 5, 5 });

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

		/**********
		 * Here I split an IterableInterval along X direction at the median
		 * X-coordinate values
		 **********/
		int n = list.numDimensions();
		PointSampleList<FloatType> LeftTreeX = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> RightTreeX = new PointSampleList<FloatType>(n);

		PointSampleList<FloatType> LeftsubTreeLeftX = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> LeftsubTreeRightX = new PointSampleList<FloatType>(n);

		PointSampleList<FloatType> RightsubTreeLeftX = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> RightsubTreeRightX = new PointSampleList<FloatType>(n);

		Pair<PointSampleList<FloatType>, PointSampleList<FloatType>> LefttreePairX = new ValuePair<PointSampleList<FloatType>, PointSampleList<FloatType>>(
				LeftsubTreeLeftX, LeftsubTreeRightX);
		Pair<PointSampleList<FloatType>, PointSampleList<FloatType>> RighttreePairX = new ValuePair<PointSampleList<FloatType>, PointSampleList<FloatType>>(
				RightsubTreeLeftX, RightsubTreeRightX);

		XcoordinatesSort = sortedCoordinates(list, 0);

		MedianLeftX = medianValueLeft(XcoordinatesSort, 0, XcoordinatesSort.size() - 1, 0);

		MedianRightX = medianValueRight(XcoordinatesSort, 0, XcoordinatesSort.size() - 1, 0);

		LeftTreeX = getLeftTree(list, MedianLeftX, 0);
		RightTreeX = getRightTree(list, MedianRightX, 0);

		for (int medianindex = 1; medianindex < MedianLeftX.size(); ++medianindex)
			LefttreePairX = getsubTrees(LeftTreeX, MedianLeftX, medianindex, 0);

		for (int medianindex = 1; medianindex < MedianRightX.size(); ++medianindex)
			RighttreePairX = getsubTrees(RightTreeX, MedianRightX, medianindex, 0);

		/*****
		 * The primary partition (along X direction) is stored in LeftTreeX and
		 * RightTreeX and partitioned space after that (also in X-direction) are
		 * stored in LefttreePairX and RighttreePairX
		 *******/

		/**********
		 * Here I repeat the above process along Y direction taking the original
		 * list for partitioning the space along Y direction
		 **********/

		PointSampleList<FloatType> LeftTreeY = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> RightTreeY = new PointSampleList<FloatType>(n);

		PointSampleList<FloatType> LeftsubTreeLeftY = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> LeftsubTreeRightY = new PointSampleList<FloatType>(n);

		PointSampleList<FloatType> RightsubTreeLeftY = new PointSampleList<FloatType>(n);
		PointSampleList<FloatType> RightsubTreeRightY = new PointSampleList<FloatType>(n);

		Pair<PointSampleList<FloatType>, PointSampleList<FloatType>> LefttreePairY = new ValuePair<PointSampleList<FloatType>, PointSampleList<FloatType>>(
				LeftsubTreeLeftY, LeftsubTreeRightY);
		Pair<PointSampleList<FloatType>, PointSampleList<FloatType>> RighttreePairY = new ValuePair<PointSampleList<FloatType>, PointSampleList<FloatType>>(
				RightsubTreeLeftY, RightsubTreeRightY);

		YcoordinatesSort = sortedCoordinates(list, 1);

		MedianLeftY = medianValueLeft(YcoordinatesSort, 0, YcoordinatesSort.size() - 1, 1);

		MedianRightY = medianValueRight(YcoordinatesSort, 0, YcoordinatesSort.size() - 1, 1);

		LeftTreeY = getLeftTree(list, MedianLeftY, 1);
		RightTreeY = getRightTree(list, MedianRightY, 1);

		for (int medianindex = 1; medianindex < MedianLeftY.size(); ++medianindex)
			LefttreePairY = getsubTrees(LeftTreeY, MedianLeftY, medianindex, 1);

		for (int medianindex = 1; medianindex < MedianRightY.size(); ++medianindex)
			RighttreePairY = getsubTrees(RightTreeY, MedianRightY, medianindex, 1);

		/*****
		 * The primary partition (along Y direction) is stored in LeftTreeY and
		 * RightTreeY and partitioned space after that (also in Y-direction) are
		 * stored in LefttreePairY and RighttreePairY
		 *******/
		double[] testlocationX;

		double[] testpoint = { 4.8, 0.4 };

		testlocationX = searchTree(MedianLeftX, MedianRightX, testpoint[0]);

		for (int i = 0; i < testlocationX.length; ++i) {
			System.out.println(testlocationX[i]);
		}
		
		double[] testlocationY;

		

		testlocationY = searchTree(MedianLeftY, MedianRightY, testpoint[1]);

		for (int i = 0; i < testlocationY.length; ++i) {
			System.out.println(testlocationY[i]);
		}
		
		
		
		System.out.println(MedianLeftX);
		System.out.println(MedianRightX);
		// System.out.println(LefttreePairX.getB().size());
		Cursor<FloatType> testtwo = LefttreePairX.getA().cursor();

		while (testtwo.hasNext()) {
			testtwo.fwd();
			Point newpoint = new Point(img.numDimensions());

			newpoint.setPosition(testtwo);

			// System.out.println("Set of x co-ordinates sorted List : " +
			// newpoint.getDoublePosition(0));
			// System.out.println("Set of y co-ordinates sorted List : " +
			// newpoint.getDoublePosition(1));
			// System.out.println("LeftTreeleft : " + testtwo.get());

		}

		Cursor<FloatType> testthree = LefttreePairX.getB().cursor();

		while (testthree.hasNext()) {
			testthree.fwd();
			Point newpointsec = new Point(img.numDimensions());

			newpointsec.setPosition(testthree);

			// System.out.println("Set of x co-ordinates sorted List : " +
			// newpointsec.getDoublePosition(0));
			// System.out.println("Set of y co-ordinates sorted List : " +
			// newpointsec.getDoublePosition(1));
			// System.out.println("LeftTreeRight : " + testthree.get());

		}

	}

}
