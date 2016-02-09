package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import net.imglib2.PointSampleList;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.kdtree.KDTreeNodeIterable;
import net.imglib2.algorithm.kdtree.SplitHyperPlaneKDTree;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTreeNode;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.TwoDtree.Distance;
import varun.TwoDtree.EucledianDistance;

public class Kdbranch {

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

	public static <T extends RealType<T>> Double medianElement(ArrayList<Long> sortedcoordinateList,
			int[] medianIndex) {

		double medianElement = 0.0;

		medianElement = 0.5 * (sortedcoordinateList.get(medianIndex[0]) + sortedcoordinateList.get(medianIndex[1]));

		return medianElement;

	}

	public static <T extends RealType<T>> int[] medianIndex(ArrayList<Long> sortedcoordinateList, int startindex,
			int lastindex, int direction) {

		// Size of a list is lastindex-startindex+1.

		int[] medianindex = new int[2];
		int medianIndexA, medianIndexB;

		if ((lastindex - startindex + 1) % 2 == 1) {
			medianIndexA = startindex + (lastindex - startindex + 1 - (lastindex - startindex + 1) % 2) / 2;
			medianIndexB = medianIndexA;
		}

		else {

			medianIndexA = startindex + (lastindex - startindex + 1) / 2 - 1;
			medianIndexB = medianIndexA + 1;
		}

		medianindex[0] = medianIndexA;

		medianindex[1] = medianIndexB;

		return medianindex;

	}

	private static class Node<T> {

		protected final double medianValue;

		protected final int direction;

		public final PointSampleList<T> LeftTree;

		public final PointSampleList<T> RightTree;

		public Node(final double medianValue, final int direction, final PointSampleList<T> LeftTree,
				final PointSampleList<T> RightTree) {

			this.medianValue = medianValue;
			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;
		}

	}

	public static <T extends RealType<T>> Node<T> getTree(PointSampleList<T> list, ArrayList<Long> sortedcoordinateList,
			int startindex, int lastindex, int direction) {
	
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

			

			if (lastindex > startindex) {

				int[] medianIndexA = new int[2];

				medianIndexA = medianIndex(sortedcoordinateList, startindex, lastindex, direction);

				double pivotElement;

				pivotElement = medianElement(sortedcoordinateList, medianIndexA);

				final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
				final PointSampleList<T> RightTree = new PointSampleList<T>(n);

				final Cursor<T> listCursor = list.localizingCursor();

				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point cord = new Point(listCursor);

					if (listCursor.getDoublePosition(direction) < pivotElement) {

						LeftTree.add(cord, listCursor.get().copy());

					} else {

						RightTree.add(cord, listCursor.get().copy());

					}

				}

			
						 
						
						 
							
							
					 
				
				
		
				
				

			
				

				return new Node<T>(pivotElement, direction, LeftTree, RightTree);

			}
			
			
			

			else
				return null;

		}

	}
	
	  public static <T extends RealType<T>> ArrayList<Node<T>> getsubTrees(Node<T> Root,
	  ArrayList<Long> sortedcoordinateList, int startindex, int lastindex, int
	  direction) {
	 
	  
	 
	  double pivotElement= Root.medianValue;
	  
	  PointSampleList<T> childA = Root.LeftTree; 
	  PointSampleList<T> childB =Root.RightTree; 
	  
	
			// Medians for left part of the tree
		 
	  ArrayList<Node<T>> allnodes = new ArrayList<Node<T>>();
		
		
		 
		if (lastindex - startindex + 1 <= 2)
			return null;
		
		else
			
	{
			
			
			  int[] medianIndexA = new int[2];
			  int[] medianIndexleftA = new int[2];
			  double medianElement=0.0;
			  double firstElement;
				
				firstElement = sortedcoordinateList.get(startindex);
			  
				medianIndexA = medianIndex(sortedcoordinateList, startindex, lastindex, direction);

				 Node<T> newnodeLeft;// = new Node<T>(pivotElement, direction, childA, childB);
				 
				 int  initialindex=medianIndexA[0] - 1;
			
			 for(int index = lastindex; index>startindex;--index){
				 
				 
				 
				 
				 
				
				 
				if (initialindex - startindex + 1 <= 2)
					break;
				
				
				 medianIndexleftA=medianIndex(sortedcoordinateList, startindex, initialindex, direction);
				
				medianElement=medianElement(sortedcoordinateList, medianIndexleftA);
				
				 if (medianElement <= firstElement)
						break;
				
				
				 newnodeLeft = getTree(childA, sortedcoordinateList, startindex, initialindex, direction);
				 
				 
				 initialindex=medianIndexleftA[0] - 1;
				 
				 
allnodes.add(newnodeLeft);				 
				 System.out.println(allnodes.get(0).LeftTree);
				 
					
					
			 }
			
			

						 
			return allnodes;
						
						
	
		
	}


	//final Node<T>	newnodeRight = getTree(childB, sortedcoordinateList, medianIndexA[0], lastindex, direction);

	 
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

	public static <T extends RealType<T>> void createBitimage(Img<T> img, Img<BitType> imgout, T ThresholdValue) {

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

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));

		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());

		FloatType val = new FloatType(100);

		// createBitimage(img, imgout, val);

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		ArrayList<Long> XcoordinatesSort = new ArrayList<Long>((int) list.dimension(0));
		ArrayList<Long> YcoordinatesSort = new ArrayList<Long>((int) list.dimension(1));

		ArrayList<Double> MedianElements = new ArrayList<Double>();

		// Make a list by setting an appropriate
		// interval on the image.

		IterableInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 10, 10 });

		final Cursor<FloatType> first = view.cursor();

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

		Node<FloatType>  testnode;
ArrayList<Node<FloatType>> allnodes = new ArrayList<Node<FloatType>>() ;
		
		
		
		int lastindex = (int) XcoordinatesSort.size() - 1;
		int startindex = 0;
		
		
		
		testnode = getTree(list, XcoordinatesSort, startindex, lastindex, 0);

		
		allnodes=getsubTrees(testnode, XcoordinatesSort, startindex, lastindex, 0);
		
		

	//	System.out.println(rootnode.LeftTree.size());
		// System.out.println(testnode.LeftTree.size());

		 Cursor<FloatType> test= allnodes.get(1).LeftTree.cursor();

		 while(test.copyCursor().hasNext()){
		 test.fwd();

		 System.out.println("allNodes index 1 lefttree : " +
		 test.get());
		 }

		/*****
		 * The primary partition (along X direction) is stored in LeftTreeX and
		 * RightTreeX and partitioned space after that (also in X-direction) are
		 * stored in LefttreePairX and RighttreePairX
		 *******/

		/**********
		 * Here I repeat the above process along Y direction taking the original
		 * list for partitioning the space along Y direction
		 **********/

		PointSampleList<BitType> LeftTreeY = new PointSampleList<BitType>(n);
		PointSampleList<BitType> RightTreeY = new PointSampleList<BitType>(n);

		YcoordinatesSort = sortedCoordinates(list, 1);

		/*****
		 * The primary partition (along Y direction) is stored in LeftTreeY and
		 * RightTreeY and partitioned space after that (also in Y-direction) are
		 * stored in LefttreePairY and RighttreePairY
		 *******/

		double[] testpoint = { 0.4, 0.44 };

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

	}

}
