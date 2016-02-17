package varun;



import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import com.sun.tools.internal.ws.processor.modeler.annotation.MakeSafeTypeVisitor;
import com.sun.tools.javac.util.Assert;

import net.imglib2.PointSampleList;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.kdtree.KDTreeNodeIterable;
import net.imglib2.algorithm.kdtree.SplitHyperPlaneKDTree;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTreeNode;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.RealLocalizable;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.neighborsearch.NearestNeighborSearch;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.TwoDtree.Distance;
import varun.TwoDtree.EucledianDistance;


public class MyKDtreeint {
	
	
	/********* For an input image returns a PointSampleList ********/
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

	/*********
	 * Starting the methods which sort an Arraylist of co-ordinates using
	 * Merge-Sort algorithm
	 *********/

	public static <T extends RealType<T>> void split(ArrayList<Long> coordinateList, int direction) {

		if (coordinateList.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) coordinateList.size() / 2;

			Iterator<Long> iterator = coordinateList.iterator();

			final ArrayList<Long> childA = new ArrayList<Long>((int) coordinateList.size() / 2);

			final ArrayList<Long> childB = new ArrayList<Long>(
					(int) (coordinateList.size() / 2 + coordinateList.size() % 2));

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
	/******** End of the Merge-Sort routine for Arraylist *********/

	/******* Returns the medianElement for input PointSampleList *******/

	public static <T extends RealType<T>> double getMedian(PointSampleList<T> list, int direction) {

		final Cursor<T> listCursor = list.localizingCursor();

		final ArrayList<Long> values = new ArrayList<Long>();

		while (listCursor.hasNext()) {
			listCursor.fwd();

			values.add(listCursor.getLongPosition(direction));

		}

		// Collections.sort(values);

		split(values, direction); // Since the list is sorted I only have to get
									// the value at the middle of the list

		int startindex = 0;
		int lastindex = values.size() - 1;

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

		double medianElement = 0.0;

		medianElement = 0.5 * (values.get(medianindex[0]) + values.get(medianindex[1]));

		return medianElement;

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	private static class searchNode<T> {

		private final int n;

		private final double medianValue;

		private final int direction;

		private final PointSampleList<T> searchBranch;

		public searchNode(final double medianValue, final int direction, final PointSampleList<T> searchBranch) {

			this.n = searchBranch.numDimensions();
			this.medianValue = medianValue;
			this.direction = direction;
			this.searchBranch = searchBranch;

		}

		public int getDirection() {
			return direction;
		}

		public double getMedianValue() {
			return medianValue;
		}

		public PointSampleList<T> getSearchBranch() {
			return searchBranch;
		}

		public int getNumdimensions() {
			return n;
		}

	}

	private static class Node<T> {

		public final int n;

		public final double medianValue;

		public final int direction;

		public final PointSampleList<T> LeftTree;

		public final PointSampleList<T> RightTree;

		public Node(final double medianValue, final int direction, final PointSampleList<T> LeftTree,
				final PointSampleList<T> RightTree) {
			assert LeftTree.numDimensions() == RightTree.numDimensions();
			this.n = LeftTree.numDimensions();
			this.medianValue = medianValue;
			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;
		}

		public int getnumDimensions() {
			return n;
		}

		public double getMedianValue() {
			return medianValue;

		}

		public int getDirection() {
			return direction;
		}

		public PointSampleList<T> getLeftTree() {
			return LeftTree;
		}

		public PointSampleList<T> getRightTree() {
			return RightTree;
		}

	}

	private static <T> void nodetoList(final Node<T> node, final ArrayList<Node<T>> allnodes) {
		allnodes.add(node);
	}

	private static <T> void searchNodetoList(final searchNode<T> searchnode, final ArrayList<searchNode<T>> allnodes) {
		allnodes.add(searchnode);
	}

	/******
	 * Returns a root tree, I do this to initialize an ArrayList<Node<T>> in the
	 * main program which I overwrite later to include all the subtrees (Clever
	 * or Dangerous?)
	 ******/

	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction) {

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

			pivotElement = getMedian(list, direction);

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

			Node<T> node = new Node<T>(pivotElement, direction, LeftTree, RightTree);

			return node;

		}

	}

	/*********
	 * Here I return an Arraylist of Node<T> type by a clever trick, I give in
	 * an ArrayList of Node<T> containing only the rootNode and in the course of
	 * the routine below overwrite that list with the nodes of all the subtrees.
	 * The list then contains the object Node<T> for all the subtrees including
	 * the rootTree. Is this clever or dangerous way to do it?
	 ********/
	public static <T extends RealType<T>> void getTree(PointSampleList<T> list, ArrayList<Node<T>> allnodes,
			int direction) {

		int n = list.numDimensions();
		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == n)
			direction = 0;
		if (list.dimension(direction) <= 2)
			return;

		else {

			// ArrayList<Node<T>> allnodes = new ArrayList<Node<T>>();

			double pivotElement;

			pivotElement = getMedian(list, direction);

			final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
			final PointSampleList<T> RightTree = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement) {

					LeftTree.add(cord, listCursor.get().copy());

				} else if (listCursor.getDoublePosition(direction) >= pivotElement) {

					RightTree.add(cord, listCursor.get().copy());

				}

			}

			Node<T> node = new Node<T>(pivotElement, direction, LeftTree, RightTree);

			int otherdirection = direction + 1;

			if (otherdirection == n)
				otherdirection = 0;

			getTree(LeftTree, allnodes, otherdirection);

			getTree(RightTree, allnodes, otherdirection);

			nodetoList(node, allnodes);

		}

	}

	/***********
	 * Returns the complete search path (splitNodes, direction of split and the
	 * search branches) giving all the nearest neighbours of a point, also
	 * stored are the farther neighbours of the testpoint,
	 * 
	 * index 0 of the nodeList stores the closest node and the search branch and
	 * lastindex stores the left or right side or the Root Tree where the first
	 * split happened to search for the point.
	 * 
	 * index 0 of the farnodeList stores the least farthest node to the search
	 * point and the lastindex stores the other side of the RootTree which
	 * should really be far far away from the given point.
	 ***********/

	public static <T extends RealType<T>> void closestNode(double[] testpoint, Node<T> Trees,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList) {

		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else

			otherdirection = direction + 1;

		double locationdiff = (testpoint[direction] - Trees.getMedianValue());

		final boolean leftbranchsearch = locationdiff < 0;

		final PointSampleList<T> searchBranch = leftbranchsearch ? Trees.LeftTree : Trees.RightTree;
		final PointSampleList<T> nonsearchBranch = leftbranchsearch ? Trees.RightTree : Trees.LeftTree;

		Node<T> nextnode;
		searchNode<T> searchnode;
		searchNode<T> nonsearchnode;

		if (searchBranch.dimension(otherdirection) > 2) {
			nextnode = makeNode(searchBranch, otherdirection);
			searchnode = new searchNode<T>(nextnode.medianValue, otherdirection, searchBranch);
			nonsearchnode = new searchNode<T>(nextnode.medianValue, otherdirection, nonsearchBranch);
		}

		else {

			nextnode = Trees;
			searchnode = new searchNode<T>(nextnode.medianValue, otherdirection, searchBranch);
			nonsearchnode = new searchNode<T>(nextnode.medianValue, otherdirection, nonsearchBranch);

		}

		if (nextnode != Trees) {

			closestNode(testpoint, nextnode, nodeList, farnodeList);
		}

		searchNodetoList(searchnode, nodeList);
		searchNodetoList(nonsearchnode, farnodeList);

	}

	public static <T extends RealType<T>> double volumeHypercube(PointSampleList<T> list, int dimensions) {
		double vol = list.dimension(dimensions);
		// dimensions (of space) = list.numDimensions() - 1;

		for (int d = dimensions - 1; d >= 0; --d)
			vol *= list.dimension(d);

		return vol;

	}

	public static <T extends RealType<T>> double NearestNeighbourSearch(double[] testpoint,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList, PointSampleList<T> list) {

		int n = list.numDimensions();

		int dimensions = n - 1;

		double Volume = volumeHypercube(list, dimensions);

		int depth = nodeList.size();

		double constantfactor = Math.sqrt(dimensions) * Math.pow(Volume, 1.0 / dimensions);

		double smallRadius = constantfactor / (Math.pow(2, depth / dimensions + 1));

		double bigRadius = constantfactor / (Math.pow(2, (depth - 1) / dimensions + 1));

		System.out.println(smallRadius);
		System.out.println(bigRadius);

		return smallRadius;

	}

	public static <T extends RealType<T>> double NearestNeighbourSearch(double[] testpoint,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList, final Distance dist) {

		double mindistance;
		double bestdistance = Double.MAX_VALUE;
		double secondbestdistance = Double.MAX_VALUE;

		final Cursor<T> listcursor = nodeList.get(0).getSearchBranch().localizingCursor();
		while (listcursor.hasNext()) {
			listcursor.fwd();
			mindistance = dist.getDistance(listcursor, testpoint);

			bestdistance = Math.min(mindistance, bestdistance);
		}
		for (int index = 1; index < nodeList.size() - 1; ++index) {

			final Cursor<T> cursor = nodeList.get(index).getSearchBranch().localizingCursor();

			while (cursor.hasNext()) {
				cursor.fwd();
				mindistance = dist.getDistance(cursor, testpoint);
				secondbestdistance = Math.min(mindistance, secondbestdistance);

				if (secondbestdistance < bestdistance) {
					bestdistance = secondbestdistance;
					continue;
				}

				else {

					break;
				}
			}
		}

		for (int index = 0; index < farnodeList.size() - 1; ++index) {

			final Cursor<T> cursor = farnodeList.get(index).getSearchBranch().localizingCursor();

			while (cursor.hasNext()) {
				cursor.fwd();
				mindistance = dist.getDistance(cursor, testpoint);
				secondbestdistance = Math.min(mindistance, secondbestdistance);

				if (secondbestdistance < bestdistance) {
					bestdistance = secondbestdistance;
					continue;
				}

				else {

					break;
				}
			}
		}

		
		RealCursor<T> c =farnodeList.get(0).getSearchBranch().cursor();
//		finalNode.getSearchBranch().cursor();

		 while (c.hasNext()) {

		c.fwd();
		System.out.println("Branch X: "+c.getDoublePosition(0));
System.out.println("Branch Y: "+c.getDoublePosition(1));

		 }
		
		
		return bestdistance;

	}

	/*************
	 * Implementation of analogous method for Lists called RetailAll() but done
	 * here for PointSampleLists to return a PointSampleList having only the
	 * common elements of two differently sized PointSampleLists, could be
	 * useful (not used currently)
	 *************/
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

	/********
	 * For a Pair of PointSampleLists, returns a single PointSampleList by
	 * combining the two pairs into one (not used currently)
	 ***********/

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

	/**********
	 * Starting the distance transform routine (not used currently)
	 **********/

	public interface Distance {

		double getDistance(Localizable cursor1, Localizable cursor2);

		<T extends RealType<T>> double getDistance(Localizable listcursor, double[] testpoint);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(Localizable cursor1, Localizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);
		}

		public <T extends RealType<T>> double getDistance(Localizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += (testpoint[d] - listcursor.getDoublePosition(d))
						* (testpoint[d] - listcursor.getDoublePosition(d));

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

		public <T extends RealType<T>> double getDistance(Localizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += Math.abs(testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	/************
	 * Creating a bitType image from an image of type T by doing thresholding
	 * (not used currently)
	 ************/

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

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		// Make a list by setting an appropriate
		// interval on the image.

		IterableInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 5, 5 });

		final Cursor<FloatType> first = view.cursor();

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(img.numDimensions());

			cord.setPosition(first);

			list.add(cord, first.get().copy());

		}

		int n = list.numDimensions();

		/********** Starting the KD-Tree creation *********/

		Node<FloatType> rootnode, finalnode;

		rootnode = makeNode(list, 0);

		ArrayList<Node<FloatType>> allnodes = new ArrayList<Node<FloatType>>();

		ArrayList<searchNode<FloatType>> searchnodes = new ArrayList<searchNode<FloatType>>();

		ArrayList<searchNode<FloatType>> nonsearchnodes = new ArrayList<searchNode<FloatType>>();

		// Make a KD-tree by splitting along the X direction first and then the
		// Y direction and so on until there is no more splitting possible
		getTree(list, allnodes, 0);

		/*********** Testing if the built KD tree is correct ***********/

		for (int i = 0; i < allnodes.size(); ++i) {
			// System.out.println("Median Value : "
			// +allnodes.get(i).medianValue);
			// System.out.println("Direction : " +allnodes.get(i).direction);
		}
		/******** Make a test point and search for the closest node *********/

		double[] testpoint = new double[2];

		testpoint[0] = 0.4;
		testpoint[1] = 1.2;

		closestNode(testpoint, rootnode, searchnodes, nonsearchnodes);

		double distance;
		distance = NearestNeighbourSearch(testpoint, searchnodes, nonsearchnodes, new EucledianDistance());

		System.out.println(distance);

		// for (int i = 0; i < searchnodes.size(); ++i) {
		// Cursor<FloatType> listcursor =
		// searchnodes.get(i).getSearchBranch().cursor();

		// while (listcursor.hasNext()) {

		// listcursor.fwd();
		// System.out.println(" list X cor:" + listcursor.getDoublePosition(0));
		// System.out.println(" list Y cor:" + listcursor.getDoublePosition(1));
		// System.out.println(" list: " + listcursor.get());

		// }

		// System.out.println(" Depth:" + searchnodes.size());

		// System.out.println(" Direction:" + searchnodes.get(i).direction);

		// System.out.println(" Search Branch:"
		// + searchnodes.get(i).searchBranch.size());

		// }
	}

}
