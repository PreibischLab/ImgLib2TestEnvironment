package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
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
import varun.MyKDtree.EucledianDistance;
import varun.MyKDtreeint.Distance;
import varun.MyKDtreeint.Node;

public class MyKDtree {

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

	public static PointSampleList<BitType> getvalueList(PointSampleList<BitType> list, int val) {

		int n = list.numDimensions();

		final Cursor<BitType> first = list.cursor();

		// A point sample list with coordinates declared and initialized.
		PointSampleList<BitType> parent = new PointSampleList<BitType>(n);

		while (first.hasNext()) {
			first.fwd();
			if (first.get().getInteger() == val) {
				Point cord = new Point(n);

				cord.setPosition(first);

				parent.add(cord, first.get().copy());

			}
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

			final Node<T> node = new Node<T>(pivotElement, direction, LeftTree, RightTree);

			int otherdirection = direction + 1;

			if (otherdirection == n)
				otherdirection = 0;

			getTree(LeftTree, allnodes, otherdirection);

			getTree(RightTree, allnodes, otherdirection);

			nodetoList(node, allnodes);

		}

	}

	/***********
	 * Returns the node closest to the testpoint stores as a single entry in the
	 * arrayList
	 ***********/

	public static <T extends RealType<T>> void closestNode(RealLocalizable testpoint, Node<T> Trees,
			ArrayList<Node<T>> singlenode, double Bestdistsquared) {

		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else
			otherdirection = direction + 1;

		double locationdiff = (testpoint.getDoublePosition(direction) - Trees.getMedianValue());

		double dist = locationdiff * locationdiff;
		double axisdiff = locationdiff * locationdiff;
		for (int d = 0; d < n; ++d) {
			if (d != direction)

				dist += testpoint.getDoublePosition(d) * testpoint.getDoublePosition(d);
		}

		Node<T> finalnode;

		final boolean leftbranchsearch = locationdiff < 0;

		final PointSampleList<T> searchBranch = leftbranchsearch ? Trees.LeftTree : Trees.RightTree;
		final PointSampleList<T> nonsearchBranch = leftbranchsearch ? Trees.RightTree : Trees.LeftTree;

		Node<T> nearnode, farnode;
		if ((searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1) > 2) {
			nearnode = makeNode(searchBranch, otherdirection);

			if ((searchBranch.realMax(direction) - searchBranch.realMin(direction) + 1) > 2)
				closestNode(testpoint, nearnode, singlenode, Bestdistsquared);

		}

		if (axisdiff <= Bestdistsquared
				&& (nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1) > 2) {
			farnode = makeNode(nonsearchBranch, otherdirection);

			if ((nonsearchBranch.realMax(direction) - nonsearchBranch.realMin(direction) + 1) > 2)
				closestNode(testpoint, farnode, singlenode, Bestdistsquared);
		}

		if (dist < Bestdistsquared) {
			Bestdistsquared = dist;

			finalnode = Trees;

			nodetoList(finalnode, singlenode);

		}

	}

	/********
	 * For a Node<T>, returns a single PointSampleList by combining the Left and
	 * Right Tree pairs into one (not used currently)
	 ***********/

	public static <T extends RealType<T>> PointSampleList<T> combineTrees(Node<T> Tree) {

		assert Tree.LeftTree.numDimensions() == Tree.RightTree.numDimensions();
		int n = Tree.LeftTree.numDimensions();

		PointSampleList<T> singleTree = new PointSampleList<T>(n);
		Cursor<T> treecursorA = Tree.LeftTree.cursor();
		Cursor<T> treecursorB = Tree.RightTree.cursor();

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

	public static <T extends RealType<T>> void testNeighbours(PointSampleList<BitType> list, final Distance dist)
			throws FileNotFoundException {

		PrintStream out = new PrintStream(new FileOutputStream("output.txt"));
		System.setOut(out);

		Node<BitType> rootnode;

		rootnode = makeNode(list, 0);
		final Cursor<BitType> listcursor = list.cursor();
		double Bestdistsquared = Double.MAX_VALUE;
		while (listcursor.hasNext()) {
			listcursor.fwd();
			ArrayList<Node<BitType>> singlenode = new ArrayList<Node<BitType>>();
			closestNode(listcursor, rootnode, singlenode, Bestdistsquared);
			PointSampleList<BitType> singletree = combineTrees(singlenode.get(0));
			Cursor<BitType> singlecursor = singletree.cursor();

			while (singlecursor.hasNext()) {

				singlecursor.fwd();

				System.out.println("Test Point X: " + listcursor.getDoublePosition(0));
				System.out.println("Test Point Y: " + listcursor.getDoublePosition(1));

				System.out.println("Neighbour Points X: " + singlecursor.getDoublePosition(0));
				System.out.println("Neighbour Points Y: " + singlecursor.getDoublePosition(1));

			}

		}

	}

	/**********
	 * Starting the distance transform routine
	 **********/

	public static <T extends RealType<T>> void distanceTransform(PointSampleList<BitType> totallist,
			PointSampleList<BitType> list, PointSampleList<BitType> listzerosorones, RandomAccessibleInterval<T> imgout,
			final Distance dist) {

		Node<BitType> rootnode;

		rootnode = makeNode(list, 0);

		double mindistance = 0;
		double Bestdistsquared = Double.MAX_VALUE;

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		while (zerooronelistcursor.hasNext()) {
			zerooronelistcursor.fwd();
			outbound.setPosition(zerooronelistcursor);
			ArrayList<Node<BitType>> singlenode = new ArrayList<Node<BitType>>();
			closestNode(zerooronelistcursor, rootnode, singlenode, Bestdistsquared);
			PointSampleList<BitType> singletree = combineTrees(singlenode.get(0));

			Cursor<BitType> singlecursor = singletree.cursor();

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				mindistance = dist.getDistance(zerooronelistcursor, singlecursor);
				outbound.get().setReal(mindistance);
				System.out.println(mindistance);
			}

		}
		
		final Cursor<BitType> listcursor = list.localizingCursor();
while(listcursor.hasNext()){
	listcursor.fwd();
	outbound.setPosition(listcursor);
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
				// if (cursor1.getDoublePosition(d) !=
				// cursor2.getDoublePosition(d))
				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);
				// else
				// distance = Double.MAX_VALUE;
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
				// if (cursor1.getDoublePosition(d) !=
				// cursor2.getDoublePosition(d))
				distance += Math.abs(cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d));
				// else
				// distance = Double.MAX_VALUE;
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

	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());

		FloatType val = new FloatType(200);

		createBitimage(img, bitimg, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		RandomAccessibleInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 50, 50 });

		list = getList(bitimg);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(bitimg.numDimensions());

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(bitimg.numDimensions());

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		 distanceTransform(list,listonlyones,listonlyzeros, imgout, new
		 EucledianDistance());

		// testNeighbours(list, new EucledianDistance()); // Writes nearest
		// neighbours in a file

		ImageJFunctions.show(imgout).setTitle("KD-Tree output");

	}
}