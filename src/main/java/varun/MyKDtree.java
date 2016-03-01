package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import com.sun.tools.internal.ws.resources.GeneratorMessages;
import com.sun.tools.javac.api.Formattable.LocalizedString;

import ij.ImageJ;
import net.imglib2.PointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTree.KDTreeCursor;
import net.imglib2.KDTreeNode;
import net.imglib2.Localizable;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPointSampleList;
import net.imglib2.Sampler;
import net.imglib2.algorithm.kdtree.KDTreeNodeIterable;
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
	public static <T extends RealType<T>> PointSampleList<T> getList(IterableInterval<T> img) {

		int n = img.numDimensions();

		final Cursor<T> first = img.cursor();

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

	public static class Node<T> {

		private final int n;

		private final double[] nodePoint;

		private final int direction;

		private final PointSampleList<T> LeftTree;

		private final PointSampleList<T> RightTree;

		private final ArrayList<Point> newleftlist;

		private final ArrayList<Point> newrightlist;

		public Node(double[] nodePoint, final int direction, final PointSampleList<T> LeftTree,
				final PointSampleList<T> RightTree, final ArrayList<Point> newleftlist, ArrayList<Point> newrightlist) {
			assert LeftTree.numDimensions() == RightTree.numDimensions();
			this.n = LeftTree.numDimensions();

			this.nodePoint = nodePoint;

			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;

			this.newleftlist = newleftlist;

			this.newrightlist = newrightlist;

		}

		public int getnumDimensions() {
			return n;
		}

		public double getMedianValue() {
			return nodePoint[direction];

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

		public void localize(double[] point) {

			for (int d = 0; d < n; ++d)
				point[d] = (double) nodePoint[d];
		}

		public ArrayList<Point> getNewleftXlist() {
			return newleftlist;
		}

		public ArrayList<Point> getNewrightXlist() {
			return newrightlist;
		}

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	public static class searchNode<T> {

		private final int n;

		public final double[] Position;

		protected final PointSampleList<T> list;

		protected final ArrayList<Point> Xlist;

		protected final ArrayList<Point> Ylist;

		protected Node<T> finalnode;

		protected ArrayList<Node<T>> allnodes = new ArrayList<Node<T>>();

		protected double Bestdistsquared;
		protected double Bestaxisdiffsquared;

		public searchNode(final PointSampleList<T> list, final ArrayList<Point> Xlist, final ArrayList<Point> Ylist) {

			n = list.numDimensions();
			Position = new double[n];
			this.Xlist = Xlist;
			this.Ylist = Ylist;
			this.list = list;

		}

		public int getNumdimensions() {
			return n;
		}

		public void search(final RealLocalizable cursor, final int direction) throws FileNotFoundException {
			cursor.localize(Position);
			Bestdistsquared = Double.MAX_VALUE;
			Bestaxisdiffsquared = Double.MAX_VALUE;
			closestNode(list, Xlist, Ylist, direction, Bestdistsquared);
		}

		public double getBestdist() {
			return Bestdistsquared;
		}

		public double getBestaxisdiffdist() {
			return Bestaxisdiffsquared;
		}

		public Node<T> getfinalnode() {
			return finalnode;
		}

		public ArrayList<Node<T>> getallnodes() {
			return allnodes;
		}

		private void closestNode(final PointSampleList<T> list, final ArrayList<Point> Xlist,
				final ArrayList<Point> Ylist, final int direction, double olddistance) throws FileNotFoundException {

			if (list.dimension(direction) <= 2) {

				return;
			}

			else {

				final boolean directionchoice = direction == n - 1;
				final int otherdirection = directionchoice ? 0 : direction + 1;

				final Node<T> currentBest = makeNode(list, Xlist, Ylist, direction);

				double dist = 0;

				for (int d = 0; d < n; ++d) {

					dist += Math.pow((Position[d] - currentBest.nodePoint[d]), 2);
				}

				final double locationdiff = Position[currentBest.direction]
						- currentBest.nodePoint[currentBest.direction];

				final double axisdiff = locationdiff * locationdiff;

				final boolean leftbranchsearch = locationdiff < 0;

				if (dist <= Bestdistsquared && axisdiff <= Bestaxisdiffsquared) {

					Bestdistsquared = dist;
					Bestaxisdiffsquared = axisdiff;
					finalnode = currentBest;

				}

				/*
				 * 
				 * if ( dist <=Bestdistsquared ) {
				 * 
				 * Bestdistsquared = dist;
				 * 
				 * finalnode = currentBest;
				 * 
				 * }
				 */

				PointSampleList<T> searchBranch = leftbranchsearch ? currentBest.LeftTree : currentBest.RightTree;
				final ArrayList<Point> sortedlist = leftbranchsearch ? currentBest.newleftlist
						: currentBest.newrightlist;

				final PointSampleList<T> nonsearchBranch = leftbranchsearch ? currentBest.RightTree
						: currentBest.LeftTree;
				final ArrayList<Point> sortedfarlist = leftbranchsearch ? currentBest.newrightlist
						: currentBest.newleftlist;

				final ArrayList<Point> newXlist, newYlist, newnonsXlist, newnonsYlist;

				final boolean nodedirectionchoice = currentBest.direction == n - 1;

				newXlist = nodedirectionchoice ? Xlist : sortedlist;

				newYlist = nodedirectionchoice ? sortedlist : Ylist;

				newnonsXlist = nodedirectionchoice ? Xlist : sortedfarlist;

				newnonsYlist = nodedirectionchoice ? sortedfarlist : Ylist;

				if (dist <= olddistance) {

					closestNode(searchBranch, newXlist, newYlist, otherdirection, dist);

					if (axisdiff <= Bestdistsquared)
						closestNode(nonsearchBranch, newnonsXlist, newnonsYlist, otherdirection, dist);

				}

			}

		}

		public double[] getMedian(ArrayList<Point> Xlist, ArrayList<Point> Ylist, int direction) {

			final boolean directionchoice = direction == n - 1;
			final int otherdirection = directionchoice ? 0 : direction + 1;

			final ArrayList<Point> cordsort = directionchoice ? Ylist : Xlist;

			final ArrayList<Point> anticordsort = directionchoice ? Xlist : Ylist;

			final double[] medianPoint = new double[n];

			int medianindexA = (cordsort.size()) / 2;

			medianPoint[direction] = (cordsort.get(medianindexA).getDoublePosition(direction));

			int medianindexC = (anticordsort.size()) / 2;

			medianPoint[otherdirection] = (anticordsort.get(medianindexC).getDoublePosition(otherdirection));

			return medianPoint;

		}

		public Node<T> makeNode(PointSampleList<T> sortedlist, ArrayList<Point> Xlist, ArrayList<Point> Ylist,
				int direction) {
			final boolean directionchoice = direction == n - 1;

			final ArrayList<Point> XorYlist = directionchoice ? Ylist : Xlist;

			double[] point = new double[n];

			point = getMedian(Xlist, Ylist, direction);

			final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
			final PointSampleList<T> RightTree = new PointSampleList<T>(n);

			final ArrayList<Point> newleftlist = new ArrayList<Point>();
			final ArrayList<Point> newrightlist = new ArrayList<Point>();

			final Cursor<T> listCursor = sortedlist.localizingCursor();

			final Iterator<Point> listIterator = XorYlist.listIterator();

			int index = 0;

			while (listIterator.hasNext()) {

				if (listIterator.next().getDoublePosition(direction) < point[direction])

					newleftlist.add(XorYlist.get(index));

				else
					newrightlist.add(XorYlist.get(index));

				index++;
			}

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(n);

				cord.setPosition(listCursor);

				if (listCursor.getDoublePosition(direction) < point[direction])

					LeftTree.add(cord, listCursor.get());

				else

					RightTree.add(cord, listCursor.get());

			}

			Node<T> node = new Node<T>(point, direction, LeftTree, RightTree, newleftlist, newrightlist);

			return node;

		}

	}

	public static <T extends RealType<T>> ArrayList<Point> getpointList(PointSampleList<T> sortedlist) {
		ArrayList<Point> pointlist = new ArrayList<Point>();

		Cursor<T> newlistCursor = sortedlist.localizingCursor();

		while (newlistCursor.hasNext()) {

			newlistCursor.fwd();

			Point newcord = new Point(newlistCursor);

			pointlist.add(newcord);

		}

		return pointlist;
	}

	/******
	 * Returns an arrayList of Points sorted in the given direction, useful for
	 * computing median
	 ******/

	public static <T extends RealType<T>> void sortpointList(ArrayList<Point> pointlist, int direction) {
		if (pointlist.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) pointlist.size() / 2;

			Iterator<Point> iterator = pointlist.iterator();

			final ArrayList<Point> childA = new ArrayList<Point>((int) pointlist.size() / 2);

			final ArrayList<Point> childB = new ArrayList<Point>((int) (pointlist.size() / 2 + pointlist.size() % 2));

			int index = 0;

			while (iterator.hasNext()) {
				iterator.next();

				if (index < splitIndex)
					childA.add(pointlist.get(index));

				else

					childB.add(pointlist.get(index));

				index++;

			}

			sortpointList(childA, direction);

			sortpointList(childB, direction);

			mergepointListValue(pointlist, childA, childB, direction);

			/********
			 * The part below removes the duplicate entries in the sorted array
			 ********/
			int j = 0;

			for (int i = 0; i < pointlist.size(); ++i) {

				j = i + 1;
				while (j < pointlist.size()) {

					if (pointlist.get(i).getDoublePosition(direction) == pointlist.get(j)
							.getDoublePosition(direction)) {

						pointlist.remove(j);

					}

					else {
						++j;
					}

				}

			}

		}

	}

	/// ***** Returns a sorted list *********////
	public static void mergepointListValue(ArrayList<Point> sortedlist, ArrayList<Point> listA, ArrayList<Point> listB,
			int direction) {

		int i = 0, j = 0, k = 0;

		while (i < listA.size() && j < listB.size()) {

			if (listA.get(i).getDoublePosition(direction) < (listB.get(j).getDoublePosition(direction))) {

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

	/******
	 * Remove duplicates from a sorted list, keeping only the first occurrence
	 ******/

	public static void removeDuplicates(ArrayList<Point> pointlist, int direction) {

		int j = 0;

		for (int i = 0; i < pointlist.size(); ++i) {

			j = i + 1;
			while (j < pointlist.size()) {

				if (pointlist.get(i).getDoublePosition(direction) == pointlist.get(j).getDoublePosition(direction)) {

					pointlist.remove(j);

				}

				else {
					++j;
				}

			}

		}

	}

	public static <T extends RealType<T>> PointSampleList<T> combineTrees(Node<T> Tree) {

		assert Tree.LeftTree.numDimensions() == Tree.RightTree.numDimensions();
		int n = Tree.LeftTree.numDimensions();

		PointSampleList<T> singleTree = new PointSampleList<T>(n);
		Cursor<T> treecursorA = Tree.LeftTree.cursor();
		Cursor<T> treecursorB = Tree.RightTree.cursor();

		while (treecursorA.hasNext()) {

			treecursorA.fwd();
			Point treepointA = new Point(n);

			treepointA.setPosition(treecursorA);

			singleTree.add(treepointA, treecursorA.get());

		}
		while (treecursorB.hasNext()) {

			treecursorB.fwd();
			Point treepointB = new Point(n);

			treepointB.setPosition(treecursorB);

			singleTree.add(treepointB, treecursorB.get());

		}

		return singleTree;

	}

	/****
	 * This is more of a "Java" way of getting NN by using an NN object and
	 * using getters and setters to get the Best Node point
	 ****/
	public static <T extends RealType<T>> void ConcisedistanceTransform(PointSampleList<BitType> list,
			ArrayList<Point> Xlist, ArrayList<Point> Ylist, PointSampleList<BitType> listzerosorones,
			RandomAccessibleInterval<T> imgout, final Distance dist) throws FileNotFoundException {

		// PrintStream out = new PrintStream(new
		// FileOutputStream("conKDtreemindist.txt"));
		// System.setOut(out);

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		final searchNode<BitType> Bestnode = new searchNode<BitType>(list, Xlist, Ylist);

		while (zerooronelistcursor.hasNext()) {

			zerooronelistcursor.fwd();
			double mindistance = Double.MAX_VALUE;

			outbound.setPosition(zerooronelistcursor);

			Bestnode.search(zerooronelistcursor, 0);

			PointSampleList<BitType> singletree = combineTrees(Bestnode.finalnode);

			Cursor<BitType> singlecursor = singletree.cursor();

			double distance = 0;

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerooronelistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			// System.out.println("Size of search node :" + singletree.size());

			// System.out.println(mindistance);

			outbound.get().setReal(mindistance);

		}

		final Cursor<BitType> listcursor = list.localizingCursor();
		while (listcursor.hasNext()) {
			listcursor.fwd();
			outbound.setPosition(listcursor);
			outbound.get().setReal(0);

		}

	}

	public static <T extends RealType<T>> ArrayList<Double> NNsearch(PointSampleList<BitType> list,
			ArrayList<Point> Xlist, ArrayList<Point> Ylist, PointSampleList<BitType> listzerosorones,
			final Distance dist) throws FileNotFoundException {

		// PrintStream out = new PrintStream(new
		// FileOutputStream("Pixelswithoutold.txt"));
		// System.setOut(out);

		final ArrayList<Double> distancelist = new ArrayList<Double>();

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final searchNode<BitType> Bestnode = new searchNode<BitType>(list, Xlist, Ylist);

		while (zerooronelistcursor.hasNext()) {

			zerooronelistcursor.fwd();
			double mindistance = Double.MAX_VALUE;

			Bestnode.search(zerooronelistcursor, 0);

			PointSampleList<BitType> singletree = combineTrees(Bestnode.finalnode);

			Cursor<BitType> singlecursor = singletree.cursor();

			double distance = 0;

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerooronelistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			distancelist.add(mindistance);

		}

		return distancelist;

	}

	public static <T extends RealType<T>> ArrayList<Double> BruteForce(IterableInterval<BitType> img,
			final Distance dist) throws FileNotFoundException {

		final Cursor<BitType> bound = img.cursor();
		final ArrayList<Double> distancelist = new ArrayList<Double>();

		while (bound.hasNext()) {
			bound.fwd();

			if (bound.get().getInteger() == 0) {
				final Cursor<BitType> second = img.cursor();
				double mindistance = Double.MAX_VALUE;

				double distance = 0;
				while (second.hasNext()) {
					if (second.next().getInteger() == 1) {

						distance = dist.getDistance(bound, second);

						mindistance = Math.min(mindistance, distance);

					}
				}

				distancelist.add(mindistance);

			}

		}
		return distancelist;
	}

	public static <T extends RealType<T>> void computeDistance(IterableInterval<BitType> img,
			RandomAccessibleInterval<T> imgout, final Distance dist) throws FileNotFoundException {

		final Cursor<BitType> bound = img.cursor();

		final RandomAccess<T> outbound = imgout.randomAccess();
		// PrintStream out = new PrintStream(new
		// FileOutputStream("BruteForcedist.txt"));
		// System.setOut(out);
		while (bound.hasNext()) {
			bound.fwd();
			outbound.setPosition(bound);

			if (bound.get().getInteger() == 0) {
				final Cursor<BitType> second = img.cursor();
				double mindistance = Double.MAX_VALUE;

				double distance = 0;
				while (second.hasNext()) {
					if (second.next().getInteger() == 1) {

						distance = dist.getDistance(bound, second);

						mindistance = Math.min(mindistance, distance);

					}
				}

				outbound.get().setReal(mindistance);

			}

			else {
				outbound.get().setReal(0);
			}

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

	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/DrosophilaWingsmall.tif"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgoutkd = new ArrayImgFactory<FloatType>().create(img, new FloatType());
		final Img<FloatType> imgoutbf = new ArrayImgFactory<FloatType>().create(img, new FloatType());
		int n = bitimg.numDimensions();

		FloatType val = new FloatType(180);

		// ImageJFunctions.show(img).setTitle("KD-Tree input");

		createBitimage(img, bitimg, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		IterableInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 127, 127 });

		list = getList(bitimg);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(n);

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(n);

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		ArrayList<Point> Xpointsort = new ArrayList<Point>((int) listonlyones.size());
		ArrayList<Point> Ypointsort = new ArrayList<Point>((int) listonlyones.size());

		Xpointsort = getpointList(listonlyones);

		Ypointsort = getpointList(listonlyones);

		sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
		sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate
		PrintStream out = new PrintStream(new FileOutputStream("LogsmallWing.txt"));
		System.setOut(out);
		long startTime = System.currentTimeMillis();

		/*******
		 * Concise DT is the method to obtain the distance transform via KD tree
		 * implementation
		 *******/

		ConcisedistanceTransform(listonlyones, Xpointsort, Ypointsort, listonlyzeros, imgoutkd,
				new EucledianDistance());
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("KD Tree time: " + totalTime);

		new ImageJ();

		ImageJFunctions.show(imgoutkd).setTitle("KD-Tree output");

		long startTimesec = System.currentTimeMillis();
		computeDistance(bitimg, imgoutbf, new EucledianDistance());
		long endTimesec = System.currentTimeMillis();
		long totalTimesec = endTimesec - startTimesec;
		System.out.println("Brute Force time: " + totalTimesec);

		new ImageJ();
		ImageJFunctions.show(imgoutbf).setTitle("Brute Force");
		/********
		 * Test of KD tree in finding the nearest neighbour versus finding them
		 * via Brute Force method
		 ********/

		ArrayList<Double> NNlist = new ArrayList<Double>();

		ArrayList<Double> BFlist = new ArrayList<Double>();

		NNlist = NNsearch(listonlyones, Xpointsort, Ypointsort, listonlyzeros, new EucledianDistance());

		BFlist = BruteForce(bitimg, new EucledianDistance());

		double count = 0;

		for (int index = 0; index < NNlist.size(); ++index) {

			if (NNlist.get(index) - BFlist.get(index) != 0) {
				count++;

			}
		}

		System.out.println("Number of non-matches:" + count);

		System.out.println("Unmatched out of: " + NNlist.size());

		double rate = count / NNlist.size();

		System.out.println("Sucess rate %: " + (1.0 - rate) * 100);

	}
}