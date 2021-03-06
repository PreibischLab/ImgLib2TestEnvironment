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
import net.imglib2.algorithm.kdtree.SplitHyperPlaneKDTree;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolator;
import net.imglib2.type.logic.BitType;
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

		split(values, direction);

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

			split(childA, direction);

			split(childB, direction);

			mergeListValue(coordinateList, childA, childB);

		}

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
			medianIndexA = startindex + (lastindex - startindex + 1 - (lastindex - startindex + 1) % 2) / 2;
			medianIndexB = medianIndexA;
		}

		else {

			medianIndexA = startindex + (lastindex - startindex + 1) / 2 - 1;
			medianIndexB = medianIndexA + 1;
		}
		
		
		int initialmedianIndex= medianIndexA-1;
		
		double medianElement;
		ArrayList<Double> leftmedians = new ArrayList<Double>();

		medianElement = 0.5 * (sortedcoordinateList.get(medianIndexA) + sortedcoordinateList.get(medianIndexB));
		leftmedians.add(medianElement);

		while (lastindex > startindex) {
// Medians for left part of the tree
			int startindexleft = startindex;
			int lastindexleft = medianIndexA - 1;
		
			
			
			int medianIndexleftA, medianIndexleftB;

			double firstElement;
			firstElement = sortedcoordinateList.get(startindex);

			if (lastindexleft - startindexleft + 1 <= 2)
				break;

			if ((lastindexleft - startindexleft + 1) % 2 == 1) {
				medianIndexleftA = startindexleft
						+ (lastindexleft - startindexleft + 1 - (lastindexleft - startindexleft + 1) % 2) / 2;
				medianIndexleftB = medianIndexleftA;
			}

			else {

				medianIndexleftA = startindexleft
						+ (lastindexleft - startindexleft + 1 + (lastindexleft - startindexleft + 1) % 2) / 2 - 1;
				medianIndexleftB = medianIndexleftA + 1;

			}

			medianElement = 0.5
					* (sortedcoordinateList.get(medianIndexleftA) + sortedcoordinateList.get(medianIndexleftB));

			
			medianIndexA = medianIndexleftA;
			
			if (medianElement > firstElement)
				leftmedians.add(medianElement);
			else 
				break;
			
			int startindexright = medianIndexleftA + 1;
			int lastindexright = initialmedianIndex;
			int medianIndexrightA, medianIndexrightB;

			double lastElement;
			lastElement = sortedcoordinateList.get(initialmedianIndex);
			if (lastindexright - startindexright + 1 <= 2)
				break;

			if ((lastindexright - startindexright + 1) % 2 == 1) {
				medianIndexrightA = startindexright
						+ (lastindexright - startindexright + 1 - (lastindexright - startindexright + 1) % 2) / 2;
				medianIndexrightB = medianIndexrightA;
			}

			else {
				medianIndexrightA = startindexright + (lastindexright - startindexright + 1) / 2 - 1;
				medianIndexrightB = medianIndexrightA + 1;

			}

			medianElement = 0.5
					* (sortedcoordinateList.get(medianIndexrightA) + sortedcoordinateList.get(medianIndexrightB));

			
			medianIndexB = medianIndexrightB;

			if (lastElement > medianElement)
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
			medianIndexA = startindex + (lastindex - startindex + 1 - (lastindex - startindex + 1) % 2) / 2;
			medianIndexB = medianIndexA;
		}

		else {

			medianIndexA = startindex + (lastindex - startindex + 1) / 2 - 1;
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
						+ (lastindexright - startindexright + 1 - (lastindexright - startindexright + 1) % 2) / 2;
				medianIndexrightB = medianIndexrightA;
			}

			else {
				medianIndexrightA = startindexright + (lastindexright - startindexright + 1) / 2 - 1;
				medianIndexrightB = medianIndexrightA + 1;

			}

			medianElement = 0.5
					* (sortedcoordinateList.get(medianIndexrightA) + sortedcoordinateList.get(medianIndexrightB));

			
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

	public static <T extends RealType<T>> Pair<PointSampleList<T>, PointSampleList<T>> getLeftsubTrees(
			PointSampleList<T> LeftorRightTree, ArrayList<Double> medianElements, int medianindex,
			int direction) {

				
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
			double pivotElement, nextElement, previousElement;
			final PointSampleList<T> childA = new PointSampleList<T>(n);
			final PointSampleList<T> childB = new PointSampleList<T>(n);

			final Cursor<T> listCursor = LeftorRightTree.localizingCursor();
			
			if (medianindex< medianElements.size()-1){
				
			pivotElement = medianElements.get(medianindex);
 
			previousElement= medianElements.get(medianindex-1);
			nextElement = medianElements.get(medianindex+1);// stopindex is the index of the next element on the last node split up tree.
			
			
			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point newpoint = new Point(n);
				newpoint.setPosition(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement && listCursor.getDoublePosition(direction) >= nextElement) {

					childA.add(newpoint, listCursor.get().copy());

				} else if (listCursor.getDoublePosition(direction) >= pivotElement
						&& listCursor.getDoublePosition(direction) <= previousElement)

				{

					childB.add(newpoint, listCursor.get().copy());

				}

			}
			}
			
			else{
				pivotElement= medianElements.get(medianElements.size()-1);
				previousElement=medianElements.get(medianElements.size()-2);
				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point newpoint = new Point(n);
					newpoint.setPosition(listCursor);

					if (listCursor.getDoublePosition(direction) < pivotElement ) {

						childA.add(newpoint, listCursor.get().copy());

					} else if (listCursor.getDoublePosition(direction) >= pivotElement
							&& listCursor.getDoublePosition(direction) <= previousElement)

					{

						childB.add(newpoint, listCursor.get().copy());

					}

				}
				
				
				
				
			}
			
			Pair<PointSampleList<T>, PointSampleList<T>> newpairA = new ValuePair<PointSampleList<T>, PointSampleList<T>>(
					childA, childB);

			

			

			return newpairA;

		}

	}

	public static <T extends RealType<T>> Pair<PointSampleList<T>, PointSampleList<T>> getRightsubTrees(
			PointSampleList<T> LeftorRightTree, ArrayList<Double> medianElements, int medianindex, 
			int direction) {

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
			double pivotElement, previousElement, nextElement;
			final PointSampleList<T> childA = new PointSampleList<T>(n);
			final PointSampleList<T> childB = new PointSampleList<T>(n);

			final Cursor<T> listCursor = LeftorRightTree.localizingCursor();

			if (medianindex< medianElements.size()-1){
			
			pivotElement = medianElements.get(medianindex);

			previousElement = medianElements.get(medianindex-1);
			
			nextElement = medianElements.get(medianindex+1);

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point newpoint = new Point(n);
				newpoint.setPosition(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement
						&& listCursor.getDoublePosition(direction) >= previousElement) {

					childA.add(newpoint, listCursor.get().copy());

				} else if (listCursor.getDoublePosition(direction) >= pivotElement
						&& listCursor.getDoublePosition(direction) <= nextElement)

				{

					childB.add(newpoint, listCursor.get().copy());

				}

			}
			}
			
			else{
				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point newpoint = new Point(n);
					newpoint.setPosition(listCursor);
					
					pivotElement=medianElements.get(medianElements.size()-1);

					previousElement=medianElements.get(medianElements.size()-2);
					
				if (listCursor.getDoublePosition(direction) < pivotElement
						&& listCursor.getDoublePosition(direction) >= previousElement) {

					childA.add(newpoint, listCursor.get().copy());

				} else if (listCursor.getDoublePosition(direction) >= pivotElement)

				{

					childB.add(newpoint, listCursor.get().copy());

				}
				
				}
				
				
			}
			
			
			Pair<PointSampleList<T>, PointSampleList<T>> newpairB=new ValuePair<PointSampleList<T>, PointSampleList<T>>(
					childA, childB); ;

	
		

			return newpairB;

		}

	}

	// This method returns the two closest iterable intervals from the selected
	// point
	public static <T extends RealType<T>> PointSampleList<T> searchTree(PointSampleList<T> LeftTree,
			PointSampleList<T> RightTree, ArrayList<Double> medianElementsLeft, ArrayList<Double> medianElementsRight,
			double coordinate, int direction) {

		int n = LeftTree.numDimensions();

		PointSampleList<T> TreebranchA = new PointSampleList<T>(n);
		PointSampleList<T> TreebranchB = new PointSampleList<T>(n);
		Pair<PointSampleList<T>, PointSampleList<T>> Treepair;

		PointSampleList<T> Treebranch = new PointSampleList<T>(n);

		
		if (coordinate < medianElementsLeft.get(0)) {

			// The point is on the Left Tree

			// Backward search
			
			if (coordinate <= medianElementsLeft.get(medianElementsLeft.size()-1)){
				
				Treepair = getLeftsubTrees(LeftTree, medianElementsLeft, medianElementsLeft.size()-1 , direction);

				Treebranch = combineTrees(Treepair);
			}
			
			for (int medianindex= medianElementsLeft.size()-2; medianindex>0; --medianindex){
				
				if (coordinate <= medianElementsLeft.get(medianindex) && coordinate >= medianElementsLeft.get(medianindex+1)){
					
					Treepair = getLeftsubTrees(LeftTree, medianElementsLeft, medianindex, direction);

					Treebranch = combineTrees(Treepair);
					
					
				}
				
				
				
			}
			
		}
			
			

			

		else if (coordinate >= medianElementsRight.get(0)) {
			
			if (coordinate >= medianElementsRight.get(medianElementsRight.size()-1)){
				
				Treepair = getRightsubTrees(RightTree, medianElementsRight, medianElementsRight.size()-1 , direction);
				Treebranch = combineTrees(Treepair);
				
			}
			
for (int medianindex= medianElementsRight.size()-2; medianindex>0; --medianindex){
				
				if (coordinate <= medianElementsRight.get(medianindex) && coordinate >= medianElementsRight.get(medianindex-1)){
					
					Treepair = getRightsubTrees(RightTree, medianElementsRight, medianindex, direction);

					Treebranch = combineTrees(Treepair);
					
					
				}
				
				
				
			}
			
		}
			
			
		// End of back search


		return Treebranch;
	}

	public static <T extends RealType<T>> PointSampleList<T> getNeighbourhood(PointSampleList<T> branchX,
			PointSampleList<T> branchY, int direction, int otherdirection) {

		int n = branchX.numDimensions();

		PointSampleList<T> localNeighbourhood = new PointSampleList<T>(n);

		PointSampleList<T> smalllist = new PointSampleList<T>(n);

		PointSampleList<T> biglist = new PointSampleList<T>(n);

		final Cursor<T> Xcursor = branchX.localizingCursor();
		final Cursor<T> Ycursor = branchY.localizingCursor();

		if (branchY.size() > branchX.size()) {

			while (Xcursor.hasNext()) {
				Xcursor.fwd();
				Point newpoint = new Point(n);
				newpoint.setPosition(Xcursor);
				smalllist.add(newpoint, Xcursor.get().copy());
			}
		}

		else {
			while (Ycursor.hasNext()) {
				Ycursor.fwd();
				Point newpointsec = new Point(n);
				newpointsec.setPosition(Ycursor);
				smalllist.add(newpointsec, Ycursor.get().copy());
			}

		}

		if (branchY.size() > branchX.size()) {

			while (Ycursor.hasNext()) {
				Ycursor.fwd();
				Point newpoint = new Point(n);
				newpoint.setPosition(Ycursor);
				biglist.add(newpoint, Ycursor.get().copy());
			}
		}

		else {
			while (Xcursor.hasNext()) {
				Xcursor.fwd();
				Point newpointsec = new Point(n);
				newpointsec.setPosition(Xcursor);
				biglist.add(newpointsec, Xcursor.get().copy());
			}

		}

		final Cursor<T> smallcursor = smalllist.localizingCursor();
		final Cursor<T> bigcursor = biglist.localizingCursor();
		bigcursor.fwd();

		while (bigcursor.hasNext()) {

			if (smallcursor.hasNext()) {

				smallcursor.fwd();
			}

			else {
				smallcursor.reset();
				smallcursor.fwd();
				if (bigcursor.hasNext())
					bigcursor.fwd();

			}
			if (bigcursor.getDoublePosition(direction) == smallcursor.getDoublePosition(direction)
					&& bigcursor.getDoublePosition(otherdirection) == smallcursor.getDoublePosition(otherdirection)) {

				Point newpoint = new Point(n);
				newpoint.setPosition(bigcursor);

				localNeighbourhood.add(newpoint, bigcursor.get().copy());
			}

		}

		return localNeighbourhood;

	}

	public static <T extends RealType<T>> PointSampleList<T> combineTrees(
			Pair<PointSampleList<T>, PointSampleList<T>> Treepair) {

		int n = Treepair.getA().numDimensions();

		PointSampleList<T> singleTree = new PointSampleList<T>(n);
		Cursor<T> treecursorA = Treepair.getA().cursor();
		Cursor<T> treecursorB = Treepair.getB().cursor();

		while (treecursorA.hasNext()) {

			treecursorA.fwd();
			Point treepoint = new Point(n);

			treepoint.setPosition(treecursorA);

			singleTree.add(treepoint, treecursorA.get().copy());

		}
		while (treecursorB.hasNext()) {

			treecursorB.fwd();
			Point treepoint = new Point(n);

			treepoint.setPosition(treecursorB);

			singleTree.add(treepoint, treecursorB.get().copy());

		}

		return singleTree;

	}

	public interface Distance {

		double getDistance(Localizable cursor1, Localizable cursor2);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(Localizable cursor1, Localizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);
		}

	}

	public static class MannhattanDistance implements Distance {

		public double getDistance(Localizable cursor1, Localizable cursor2) {
			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.abs(cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d));

			}

			return distance;
		}

	}

	public static <T extends RealType<T>> void createBitimage(Img<T> img, Img<BitType> imgout,
			T ThresholdValue) {

		final Cursor<T> bound = img.localizingCursor();

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

	public static  void computeDistance(PointSampleList<BitType> list,
			ArrayList<Double> MedianLeftX, ArrayList<Double> MedianRightX, ArrayList<Double> MedianLeftY,
			ArrayList<Double> MedianRightY,RandomAccessibleInterval<BitType> imgout, int direction, int otherdirection, final Distance dist) {

		PointSampleList<BitType> LeftTreeX;
		PointSampleList<BitType> RightTreeX;
		PointSampleList<BitType> LeftTreeY;
		PointSampleList<BitType> RightTreeY;
		PointSampleList<BitType> TreepairX;
		PointSampleList<BitType> TreepairY;

		LeftTreeX = getLeftTree(list, MedianLeftX, 0);
		RightTreeX = getRightTree(list, MedianRightX, 0);
		LeftTreeY = getLeftTree(list, MedianLeftY, 1);
		RightTreeY = getRightTree(list, MedianRightY, 1);

		PointSampleList<BitType> Neighbourhood;

		final Cursor<BitType> listcursor = list.localizingCursor();
		final RandomAccess<BitType> outbound = imgout.randomAccess();
		
		double distance;
		

		while (listcursor.hasNext()) {

			listcursor.fwd();
			
			outbound.setPosition(listcursor);

			if (listcursor.get().getInteger() == 0) {
				

				TreepairX = searchTree(LeftTreeX, RightTreeX, MedianLeftX, MedianRightX,listcursor.getDoublePosition(direction) , direction);
				TreepairY = searchTree(LeftTreeY, RightTreeY, MedianLeftY, MedianRightY, listcursor.getDoublePosition(otherdirection), otherdirection);
				Neighbourhood = getNeighbourhood(TreepairX, TreepairY, 0, 1);

				final Cursor<BitType> searchcursor = Neighbourhood.localizingCursor();
				double mindistance = Double.MAX_VALUE;
				while (searchcursor.hasNext()) {

					if (searchcursor.next().getInteger() == 1) {
						distance = dist.getDistance(listcursor, searchcursor);

						mindistance = Math.min(mindistance, distance);

					}
				}
				outbound.get().setReal(mindistance);
			}
			else
				outbound.get().setReal(0);
		
		}

	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));

		final Img<BitType> imgout = new ArrayImgFactory<BitType>().create(img, new BitType());

		FloatType val = new FloatType(100);

		createBitimage(img, imgout, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(img.numDimensions());

		ArrayList<Long> XcoordinatesSort = new ArrayList<Long>((int) list.dimension(0));
		ArrayList<Long> YcoordinatesSort = new ArrayList<Long>((int) list.dimension(1));

		ArrayList<Double> MedianLeftX = new ArrayList<Double>();

		ArrayList<Double> MedianRightX = new ArrayList<Double>();

		ArrayList<Double> MedianLeftY = new ArrayList<Double>();

		ArrayList<Double> MedianRightY = new ArrayList<Double>();

		// Make a list by setting an appropriate
		// interval on the image.

		IterableInterval<BitType> view = Views.interval(imgout, new long[] { 0, 0 }, new long[] { 20, 20 });
		

		final Cursor<BitType> first = imgout.cursor();

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(img.numDimensions());

			cord.setPosition(first);

			list.add(cord, first.get().copy());

		}

		/**********
		 * Here I split an IterableInterval along X direction at the median
		 * X-coordinate values
		 **********/
		int n = list.numDimensions();
		PointSampleList<BitType> LeftTreeX = new PointSampleList<BitType>(n);
		PointSampleList<BitType> RightTreeX = new PointSampleList<BitType>(n);

		XcoordinatesSort = sortedCoordinates(list, 0);
		
		System.out.println(XcoordinatesSort.size());

		MedianLeftX = medianValueLeft(XcoordinatesSort, 0, XcoordinatesSort.size() - 1, 0);

		
		System.out.println(MedianLeftX);
		/*
		MedianRightX = medianValueRight(XcoordinatesSort, 0, XcoordinatesSort.size() - 1, 0);

		LeftTreeX = getLeftTree(list, MedianLeftX, 0);
		RightTreeX = getRightTree(list, MedianRightX, 0);

		System.out.println(MedianLeftX);

		/*****
		 * The primary partition (along X direction) is stored in LeftTreeX and
		 * RightTreeX and partitioned space after that (also in X-direction) are
		 * stored in LefttreePairX and RighttreePairX
		 *******/

		/**********
		 * Here I repeat the above process along Y direction taking the original
		 * list for partitioning the space along Y direction
		 **********/
  /*
		PointSampleList<BitType> LeftTreeY = new PointSampleList<BitType>(n);
		PointSampleList<BitType> RightTreeY = new PointSampleList<BitType>(n);

	//	YcoordinatesSort = sortedCoordinates(list, 1);

	//	MedianLeftY = medianValueLeft(YcoordinatesSort, 0, YcoordinatesSort.size() - 1, 1);

	//	MedianRightY = medianValueRight(YcoordinatesSort, 0, YcoordinatesSort.size() - 1, 1);

	//	LeftTreeY = getLeftTree(list, MedianLeftY, 1);
	//	RightTreeY = getRightTree(list, MedianRightY, 1);

		*/
		/*****
		 * The primary partition (along Y direction) is stored in LeftTreeY and
		 * RightTreeY and partitioned space after that (also in Y-direction) are
		 * stored in LefttreePairY and RighttreePairY
		 *******/

	//	double[] testpoint = { 0.4, 0.44 };

		/*****
		 * Do this if you want to return only the immediate neighborhood of the
		 * point: PointSampleList<FloatType> TreebranchX; PointSampleList
		 * <FloatType> TreebranchY;
		 * 
		 * TreebranchX = searchTree(LeftTreeX, RightTreeX, MedianLeftX,
		 * MedianRightX, testpoint[0], 0);
		 * 
		 * TreebranchY = searchTree(LeftTreeY, RightTreeY, MedianLeftY,
		 * MedianRightY, testpoint[1], 1);
		 *****/

		/**********
		 * Do this if you want to return also one level up of the neighborhood
		 * of the point
		 **************/

	//	PointSampleList<BitType> TreepairX;
	//	PointSampleList<BitType> TreepairY;

	//	TreepairX = searchTree(LeftTreeX, RightTreeX, MedianLeftX, MedianRightX, testpoint[0], 0);
	//	TreepairY = searchTree(LeftTreeY, RightTreeY, MedianLeftY, MedianRightY, testpoint[1], 1);

	//	PointSampleList<BitType> Neighbourhood = new PointSampleList<BitType>(n);

	//	Neighbourhood = getNeighbourhood(TreepairX, TreepairY, 0, 1);

		//Cursor<BitType> testfour = Neighbourhood.cursor();

		
//final Img<BitType> bitimgout = new ArrayImgFactory<BitType>().create(imgout, new BitType());
		
 //computeDistance(list,
//		MedianLeftX, MedianRightX, MedianLeftY,
//		MedianRightY,bitimgout, 
//	0, 1, new EucledianDistance());


	//	Cursor<BitType> test= bitimgout.cursor();
		
	//	while (test.hasNext()) {
	//		test.fwd();
		
		//	Point newpointthird = new Point(img.numDimensions());

		//	newpointthird.setPosition(test);

		//	System.out.println("Set of x co-ordinates sorted List : " + newpointthird.getDoublePosition(0));
		//	System.out.println("Set of y co-ordinates sorted List : " + newpointthird.getDoublePosition(1));
		//	 System.out.println("Neighbourhood having the point : " +
		//	 test.get());

	//	}

		
	//	ImageJFunctions.show(bitimgout);
		
		
		
	}

}
