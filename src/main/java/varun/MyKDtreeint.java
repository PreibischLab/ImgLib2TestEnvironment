package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import net.imglib2.PointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPointSampleList;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.MyKDtree.Distance;
import varun.MyKDtree.EucledianDistance;
import varun.MyKDtree.Node;
import varun.MyKDtree.searchNode;

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

	public static <T extends RealType<T>> void split(ArrayList<Double> coordinateList, int direction) {

		if (coordinateList.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) coordinateList.size() / 2;

			Iterator<Double> iterator = coordinateList.iterator();

			final ArrayList<Double> childA = new ArrayList<Double>((int) coordinateList.size() / 2);

			final ArrayList<Double> childB = new ArrayList<Double>(
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
	public static <T extends RealType<T>> void mergeListValue(ArrayList<Double> sortedlist, ArrayList<Double> listA,
			ArrayList<Double> listB) {

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

		final ArrayList<Double> values = new ArrayList<Double>();

		while (listCursor.hasNext()) {
			listCursor.fwd();

			values.add(listCursor.getDoublePosition(direction));

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
//System.out.println(medianIndexA);
		medianElement = 0.5 * (values.get(medianindex[0]) + values.get(medianindex[1]));

		return medianElement;

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	public static class searchNode<T> {

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

	public static class Node<T> {

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
		if ((list.realMax(direction) - list.realMin(direction) + 1) <= 2)
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

				if (listCursor.getDoublePosition(direction) < pivotElement)

					LeftTree.add(cord, listCursor.get().copy());

				else

					RightTree.add(cord, listCursor.get().copy());

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
		if ((list.realMax(direction) - list.realMin(direction) + 1) <= 2)
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

				if (listCursor.getDoublePosition(direction) < pivotElement)

					LeftTree.add(cord, listCursor.get().copy());

				else

					RightTree.add(cord, listCursor.get().copy());

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

	public static <T extends RealType<T>> void closestNode(double[] testpoint, ArrayList<Node<T>> Trees,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList) {

		for (int index = Trees.size() - 1; index >= 0; --index) {
			int direction = Trees.get(index).direction;

			double locationdiff = (testpoint[direction] - Trees.get(index).getMedianValue());

			final boolean rightbranchsearch = locationdiff >= 0;

			final PointSampleList<T> searchBranch = rightbranchsearch ? Trees.get(index).RightTree
					: Trees.get(index).LeftTree;
			final PointSampleList<T> nonsearchBranch = rightbranchsearch ? Trees.get(index).LeftTree
					: Trees.get(index).RightTree;

			final searchNode<T> searchnode = new searchNode<T>(Trees.get(index).medianValue, direction, searchBranch);
			final searchNode<T> nonsearchnode = new searchNode<T>(Trees.get(index).medianValue, direction,
					nonsearchBranch);

			searchNodetoList(searchnode, nodeList);
			searchNodetoList(nonsearchnode, farnodeList);

		}

	}

	public static <T extends RealType<T>> void closestNode(RealLocalizable testpoint, Node<T> searchnode,
			Node<T> farnode) {

		int n = searchnode.getnumDimensions();
		int otherdirection;

		int direction = searchnode.direction;

		if (direction == n - 1)
			otherdirection = 0;

		else
			otherdirection = direction + 1;

		double locationdiff = (testpoint.getDoublePosition(direction) - searchnode.getMedianValue());

		final boolean rightbranchsearch = locationdiff >= 0;

		final PointSampleList<T> searchBranch = rightbranchsearch ? searchnode.RightTree : searchnode.LeftTree;
		final PointSampleList<T> nonsearchBranch = rightbranchsearch ? searchnode.LeftTree : searchnode.RightTree;

		if (searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1 > 2) {

			searchnode = makeNode(searchBranch, otherdirection);

		}
		if (nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1 > 2) {

			farnode = makeNode(nonsearchBranch, otherdirection);

		}
		if (searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1 > 2
				|| nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1 > 2)
			closestNode(testpoint, searchnode, farnode);
		


	}

	public static <T extends RealType<T>> Pair<Double, Node<T>> NearestNeighbourSearch(RealLocalizable testpoint,
			Node<T> searchnode, Node<T> farnode, final Distance dist) {

		double mindistance;
		double bestdistance = Double.MAX_VALUE;
		double secondbestdistance = Double.MAX_VALUE;

		Node<T> finalNode = searchnode;

		final Cursor<T> listcursorleft = searchnode.getLeftTree().localizingCursor();

		while (listcursorleft.hasNext()) {
			listcursorleft.fwd();

			mindistance = dist.getDistance(listcursorleft, testpoint);

			bestdistance = Math.min(mindistance, bestdistance);

			finalNode = searchnode;

		}

		final Cursor<T> listcursorright = searchnode.getRightTree().localizingCursor();

		while (listcursorright.hasNext()) {
			listcursorright.fwd();

			mindistance = dist.getDistance(listcursorright, testpoint);

			secondbestdistance = Math.min(mindistance, bestdistance);

			if (secondbestdistance < bestdistance) {
				bestdistance = secondbestdistance;
				finalNode = searchnode;
			}

		}

		final Cursor<T> farcursorleft = farnode.getLeftTree().localizingCursor();

		while (farcursorleft.hasNext()) {
			farcursorleft.fwd();

			mindistance = dist.getDistance(farcursorleft, testpoint);
			secondbestdistance = Math.min(mindistance, bestdistance);
			if (secondbestdistance < bestdistance) {
				bestdistance = secondbestdistance;
				finalNode = searchnode;
			}

		}

		final Cursor<T> farcursorright = farnode.getRightTree().localizingCursor();

		while (farcursorright.hasNext()) {
			farcursorright.fwd();

			mindistance = dist.getDistance(farcursorright, testpoint);
			secondbestdistance = Math.min(mindistance, bestdistance);
			if (secondbestdistance < bestdistance) {
				bestdistance = secondbestdistance;
				finalNode = searchnode;
			}

		}

		Pair<Double, Node<T>> pair = new ValuePair<Double, Node<T>>(bestdistance, finalNode);

		return pair;

	}

	public static double ValueNeighbourSearch(RealLocalizable testpoint, Node<BitType> searchnode,
			Node<BitType> farnode, final Distance dist) {

		double mindistance;
		double bestdistance = Double.MAX_VALUE;
		double secondbestdistance = Double.MAX_VALUE;

		final Cursor<BitType> listcursorleft = searchnode.getLeftTree().localizingCursor();

		while (listcursorleft.hasNext()) {
			listcursorleft.fwd();
			if (listcursorleft.get().getInteger() == 1) {
				mindistance = dist.getDistance(listcursorleft, testpoint);

				bestdistance = Math.min(mindistance, bestdistance);

			}

		}

		final Cursor<BitType> listcursorright = searchnode.getRightTree().localizingCursor();

		while (listcursorright.hasNext()) {
			listcursorright.fwd();
			if (listcursorright.get().getInteger() == 1) {
				mindistance = dist.getDistance(listcursorright, testpoint);

				secondbestdistance = Math.min(mindistance, bestdistance);

				if (secondbestdistance < bestdistance)
					bestdistance = secondbestdistance;

			}

		}

		final Cursor<BitType> farcursorleft = farnode.getLeftTree().localizingCursor();

		while (farcursorleft.hasNext()) {
			farcursorleft.fwd();
			if (farcursorleft.get().getInteger() == 1) {
				mindistance = dist.getDistance(farcursorleft, testpoint);
				secondbestdistance = Math.min(mindistance, bestdistance);
				if (secondbestdistance < bestdistance)
					bestdistance = secondbestdistance;

			}

		}

		final Cursor<BitType> farcursorright = farnode.getRightTree().localizingCursor();

		while (farcursorright.hasNext()) {
			farcursorright.fwd();
			if (farcursorright.get().getInteger() == 1) {
				mindistance = dist.getDistance(farcursorright, testpoint);
				secondbestdistance = Math.min(mindistance, bestdistance);
				if (secondbestdistance < bestdistance)
					bestdistance = secondbestdistance;

			}

		}

		return bestdistance;

	}

	/**********
	 * Starting the distance transform routine
	 **********/

	public static <T extends RealType<T>> void distanceTransform(PointSampleList<BitType> list,
			RandomAccessibleInterval<T> imgout, final Distance dist) {

		int n = list.numDimensions();

		Node<BitType> rootnode, searchnode, nonsearchnode;

		ArrayList<Node<BitType>> allnodes = new ArrayList<Node<BitType>>();

		getTree(list, allnodes, 0);

		rootnode = makeNode(list, 0);

		searchnode = rootnode;
		nonsearchnode = rootnode;

		double distance;

		final Cursor<BitType> listcursor = list.cursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		while (listcursor.hasNext()) {
			listcursor.fwd();

			outbound.setPosition(listcursor);

			if (listcursor.get().getInteger() == 0) {

				closestNode(listcursor, searchnode, nonsearchnode);
				
				
				distance = ValueNeighbourSearch(listcursor, searchnode, nonsearchnode, dist);

				
				outbound.get().setReal(distance);

			}

			else

				outbound.get().setReal(0);

		}

	}

	public interface Distance {

		double getDistance(RealLocalizable cursor1, RealLocalizable cursor2);

		<T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);
		}

		public <T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += (testpoint[d] - listcursor.getDoublePosition(d))
						* (testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	public static class MannhattanDistance implements Distance {

		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {
			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.abs(cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d));

			}

			return distance;
		}

		public <T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint) {
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

	public static <T extends RealType<T>> void createBitimage(RandomAccessibleInterval<T> img, Img<BitType> imgout,
			T ThresholdValue) {

		final Cursor<T> bound = Views.iterable(img).localizingCursor();

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

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());

		FloatType val = new FloatType(200);

		createBitimage(img, bitimg, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		list = getList(bitimg);

		distanceTransform(list, imgout, new EucledianDistance());

		ImageJFunctions.show(imgout).setTitle("KD-Tree output");

	}
}