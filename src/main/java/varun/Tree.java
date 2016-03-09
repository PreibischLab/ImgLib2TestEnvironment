package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import javax.management.ImmutableDescriptor;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Tree {

	public static class Node<T> {
		private final int n;

		private final double[] nodePoint;

		private final int direction;

		private final int treeindex;

		public final PointSampleList<T> Leftsublist;

		public final PointSampleList<T> Rightsublist;

		public Node<T> leftchild;

		public Node<T> rightchild;

		public Node(double[] nodePoint, final int direction, final int treeindex, final PointSampleList<T> Leftsublist,
				final PointSampleList<T> Rightsublist, final Node<T> leftchild, final Node<T> rightchild) {

			assert Leftsublist.numDimensions() == Rightsublist.numDimensions();
			this.n = Rightsublist.numDimensions();

			this.nodePoint = nodePoint;

			this.treeindex = treeindex;
			this.direction = direction;
			this.Leftsublist = Leftsublist;
			this.Rightsublist = Rightsublist;
			this.leftchild = leftchild;
			this.rightchild = rightchild;

		}

		public Node(final Node<T> node) {
			this.n = node.n;

			this.nodePoint = node.nodePoint;
			this.direction = node.direction;
			this.Leftsublist = node.Leftsublist;
			this.Rightsublist = node.Rightsublist;
			this.leftchild = node.leftchild;
			this.rightchild = node.rightchild;
			this.treeindex = node.treeindex;

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

		public PointSampleList<T> getLeftsublist() {
			return Leftsublist;
		}

		public PointSampleList<T> getRightsublist() {
			return Rightsublist;
		}

		public void localize(double[] point) {

			for (int d = 0; d < n; ++d)
				point[d] = (double) nodePoint[d];
		}

		public int getTreeindex() {
			return treeindex;
		}

	}

	public static class searchNode<T extends RealType<T>> {

		private final int n;

		public final double[] Position;

		protected Node<T> treeRoot;
		
		protected Node<T> finalnode;

		protected double Bestdistsquared;

		public searchNode(final Node<T> treeRoot) {

			n = treeRoot.getnumDimensions();
			Position = new double[n];
			this.treeRoot = treeRoot;

		}

		public int getNumdimensions() {
			return n;
		}

		public Node<T> search(final RealLocalizable cursor) throws FileNotFoundException {

			cursor.localize(Position);

			return closestNode(treeRoot, 0);
		}

		
		public void anothersearch(final RealLocalizable cursor) throws FileNotFoundException {
			cursor.localize(Position);
			Bestdistsquared = Double.MAX_VALUE;

			closestNode(treeRoot);
		}
		
		public double getBestdist() {
			return Bestdistsquared;
		}

		public Node<T> getfinalnode() {
			return finalnode;
		}
		
		public final double sqDist(Node<T> Node) {

			double distance = 0.0;

			for (int d = 0; d < n; ++d) {

			//	 if (d== Node.direction)

				distance+= (Position[d] - Node.nodePoint[d]) * (Position[d] - Node.nodePoint[d]);

			//	 else

			//	 distance+=Position[d]*Position[d];
			}

			return distance;
		}

		private Node<T> closestNode(final Node<T> currentNode, final int treeindex) throws FileNotFoundException {

			final double locationdiff = Position[currentNode.direction] - currentNode.nodePoint[currentNode.direction];

			double axisdiff = locationdiff * locationdiff;

			final boolean leftbranchsearch = locationdiff < 0;

			final Node<T> searchBranch = leftbranchsearch ? currentNode.leftchild : currentNode.rightchild;

			final Node<T> nonsearchBranch = leftbranchsearch ? currentNode.rightchild : currentNode.leftchild;

			// Goes to the Leaf Node
			Node<T> BestNode = (searchBranch == null) ? currentNode : closestNode(searchBranch, treeindex + 1);

			final double dist = sqDist(currentNode);

			final double bestdist = sqDist(BestNode);

			if (dist <= bestdist)

				BestNode = currentNode;

			if (nonsearchBranch != null) {

				if (axisdiff <= bestdist) {

					Node<T> possibleBest = closestNode(nonsearchBranch, treeindex + 1);

					final double possiblebestdist = sqDist(possibleBest);

					if (possiblebestdist <= bestdist)
						BestNode = possibleBest;

				}
			}

			return BestNode;

		}
		
		
		
		private void closestNode(final Node<T> currentBest) throws FileNotFoundException {

			

			final double dist = sqDist(currentBest);
			
			
		
			
			final double locationdiff = Position[currentBest.direction] - currentBest.nodePoint[currentBest.direction];

			double axisdiff = locationdiff * locationdiff;

			final boolean leftbranchsearch = locationdiff < 0;

			if (dist < Bestdistsquared) {

				Bestdistsquared = dist;
				finalnode = currentBest;

			}

			final Node<T> searchBranch = leftbranchsearch ? currentBest.leftchild : currentBest.rightchild;

			final Node<T> nonsearchBranch = leftbranchsearch ? currentBest.rightchild : currentBest.leftchild;

			if (searchBranch != null)

				closestNode(searchBranch);

			if (nonsearchBranch != null && axisdiff <= Bestdistsquared)

				closestNode(nonsearchBranch);

		}
		

	}

	public static <T extends RealType<T>> double[] getMedian(ArrayList<Point> cordsort, int direction, int n) {

		final boolean directionchoice = direction == n - 1;
		final int otherdirection = directionchoice ? 0 : direction + 1;

		final double[] medianPoint = new double[n];

		int medianindexA = (cordsort.size()) / 2;
		if (cordsort.size() > 2) {
			medianPoint[direction] = cordsort.get(medianindexA).getDoublePosition(direction);
			medianPoint[otherdirection] = cordsort.get(medianindexA).getDoublePosition(otherdirection);

			//System.out.println(medianPoint[direction]);
			//System.out.println(medianPoint[otherdirection]);
			
			return medianPoint;
		}

		else
			return null;

	}

	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction, int treeindex) {
		int n = list.numDimensions();

		final boolean directionchoice = direction == n - 1;
		final int otherdirection = directionchoice ? 0 : direction + 1;

		ArrayList<Point> Xpointsort = new ArrayList<Point>();
		ArrayList<Point> Ypointsort = new ArrayList<Point>();

		Xpointsort = getpointList(list);

		Ypointsort = getpointList(list);

		//System.out.println("Xpoint:"+Xpointsort);
		//System.out.println("Ypoint:"+Ypointsort);
		
		sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
		sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate

		
		
		final ArrayList<Point> cordsort = directionchoice ? Ypointsort : Xpointsort;

		final ArrayList<Point> anticordsort = directionchoice ? Xpointsort : Ypointsort;

		
		//System.out.println("sorted points:"+cordsort);
		
		double[] point = new double[n];

		point = getMedian(cordsort, direction, n);
		
		
		

		final PointSampleList<T> Leftsublist = new PointSampleList<T>(n);
		final PointSampleList<T> Rightsublist = new PointSampleList<T>(n);

		final Cursor<T> listCursor = list.localizingCursor();
		if (point != null) {
			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

	if (listCursor.getDoublePosition(direction) < point[direction])

		Leftsublist.add(cord, listCursor.get());
	
	else
		
		Rightsublist.add(cord, listCursor.get());
	
}


				double [] pivotpoint = new double [n];
				
				Cursor<T> cursor = Rightsublist.cursor();
				cursor.fwd();

					cursor.localize(pivotpoint);

			

			final Node<T> node = new Node<T>(pivotpoint, direction, treeindex, Leftsublist, Rightsublist, null, null);

			if (node.Leftsublist.realMax(otherdirection) - node.Leftsublist.realMin(otherdirection) > 0 && node != null)

				node.leftchild = makeNode(node.Leftsublist, otherdirection, treeindex + 1);
		
			if (node.Rightsublist.realMax(otherdirection) - node.Rightsublist.realMin(otherdirection) > 0
					&& node != null)
				
				node.rightchild = makeNode(node.Rightsublist, otherdirection, treeindex + 1);
			
	
			
			return node;
		}

		else

			return null;

	}

	public static <T extends RealType<T>> ArrayList<Double> ConcisedistanceTransform(final Node<BitType> rootnode,
			PointSampleList<BitType> list, RandomAccessibleInterval<T> imgout, final Distance dist)
			throws FileNotFoundException {

		final ArrayList<Double> distancelist = new ArrayList<Double>();
		final Cursor<BitType> zerolistcursor = list.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		// This is the tree of 1's.

		final searchNode<BitType> Bestnode = new searchNode<BitType>(rootnode);

		zerolistcursor.reset();

		while (zerolistcursor.hasNext()) {

			zerolistcursor.fwd();
			double mindistance = Double.MAX_VALUE;

			outbound.setPosition(zerolistcursor);

			final Node<BitType> finalnode = Bestnode.search(zerolistcursor);

			PointSampleList<BitType> singletree = combineTrees(finalnode);

			Cursor<BitType> singlecursor = singletree.cursor();

			double distance = 0;

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerolistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			distancelist.add(mindistance);
			outbound.get().setReal(mindistance);

		}

		return distancelist;

	}

	public static <T extends RealType<T>> ArrayList<Double> TestConcisedistanceTransform(final Node<BitType> rootnode,
			PointSampleList<BitType> list, RandomAccessibleInterval<T> imgout, final Distance dist)
			throws FileNotFoundException {

		final ArrayList<Double> distancelist = new ArrayList<Double>();
		final Cursor<BitType> listcursor = list.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		// This is the tree of 1's.

		final searchNode<BitType> Bestnode = new searchNode<BitType>(rootnode);
		listcursor.reset();
		while (listcursor.hasNext()) {

			listcursor.fwd();
			double mindistance = Double.MAX_VALUE;

			outbound.setPosition(listcursor);

			if (listcursor.get().getInteger() == 0) {

				// System.out.println("Test point X:"+
				// listcursor.getDoublePosition(0));
				// System.out.println("Test point Y:"+
				// listcursor.getDoublePosition(1));

			//	final Node<BitType> finalnode = Bestnode.search(listcursor);

				Bestnode.anothersearch(listcursor);
				
				PointSampleList<BitType> singletree = combineTrees(Bestnode.finalnode);

				Cursor<BitType> singlecursor = singletree.cursor();

				double distance = 0;

				while (singlecursor.hasNext()) {
					singlecursor.fwd();

					distance = dist.getDistance(listcursor, singlecursor);
					mindistance = Math.min(distance, mindistance);

				}

				distancelist.add(mindistance);
				outbound.get().setReal(mindistance);
			}

			else
				outbound.get().setReal(0);

		}

		return distancelist;

	}

	public static <T extends RealType<T>> PointSampleList<T> combineTrees(Node<T> Tree) {

		int n;

		if (Tree.Leftsublist != null)
			n = Tree.Leftsublist.numDimensions();
		else
			n = Tree.Rightsublist.numDimensions();

		PointSampleList<T> singleTree = new PointSampleList<T>(n);

		if (Tree.Leftsublist != null) {
			Cursor<T> treecursorA = Tree.Leftsublist.cursor();

			while (treecursorA.hasNext()) {

				treecursorA.fwd();
				Point treepointA = new Point(treecursorA);

				singleTree.add(treepointA, treecursorA.get());

			}

		}

		if (Tree.Rightsublist != null) {

			Cursor<T> treecursorB = Tree.Rightsublist.cursor();
			while (treecursorB.hasNext()) {

				treecursorB.fwd();
				Point treepointB = new Point(treecursorB);

				singleTree.add(treepointB, treecursorB.get());
				// System.out.println("X:"+treecursorB.getDoublePosition(0));
				// System.out.println("Y:"+treecursorB.getDoublePosition(1));

			}

		}

		return singleTree;

	}

	public static <T extends RealType<T>> ArrayList<Double> BruteForce(PointSampleList<BitType> listzeros,
			PointSampleList<BitType> listones, RandomAccessibleInterval<T> imgout, final Distance dist)
			throws FileNotFoundException {

		ArrayList<Double> distancelist = new ArrayList<Double>();

		final Cursor<BitType> bound = listzeros.cursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		while (bound.hasNext()) {
			bound.fwd();
			outbound.setPosition(bound);

			final Cursor<BitType> second = listones.cursor();
			double mindistance = Double.MAX_VALUE;

			double distance = 0;
			while (second.hasNext()) {
				second.fwd();

				distance = dist.getDistance(bound, second);

				mindistance = Math.min(mindistance, distance);

			}

			outbound.get().setReal(mindistance);
			distancelist.add(mindistance);

		}
		return distancelist;

	}

	public static <T extends RealType<T>> ArrayList<Double> TestBruteForce(PointSampleList<BitType> list,
			PointSampleList<BitType> listones, RandomAccessibleInterval<T> imgout, final Distance dist)
			throws FileNotFoundException {

		ArrayList<Double> distancelist = new ArrayList<Double>();

		final Cursor<BitType> bound = list.cursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		while (bound.hasNext()) {
			bound.fwd();
			outbound.setPosition(bound);
			if (bound.get().getInteger() == 0) {
				final Cursor<BitType> second = listones.cursor();
				double mindistance = Double.MAX_VALUE;

				double distance = 0;
				while (second.hasNext()) {
					second.fwd();

					distance = dist.getDistance(bound, second);

					mindistance = Math.min(mindistance, distance);

				}

				outbound.get().setReal(mindistance);
				distancelist.add(mindistance);
			} else
				outbound.get().setReal(0);

		}
		return distancelist;

	}

	public static void sortpointList(ArrayList<Point> pointlist, int direction) {
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
			 * The part below removes the duplicate entries in the sorted array (keeps the point with minimum value in other direction)
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

			if (first.next().getInteger() == val) {
				Point cord = new Point(n);

				cord.setPosition(first);

				parent.add(cord, first.get());

			}
		}
		return parent;

	}

	public interface Distance {

		double getDistance(RealLocalizable cursor1, RealLocalizable cursor2);

		double getDistance(RealLocalizable listcursor, double[] testpoint);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);

		}

		public double getDistance(RealLocalizable listcursor, double[] testpoint) {
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

		public double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += Math.abs(testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());
		final Img<FloatType> brimgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());
		int n = bitimg.numDimensions();

		ArrayList<Double> bruteforce = new ArrayList<Double>();
		ArrayList<Double> kdtree = new ArrayList<Double>();

		FloatType val = new FloatType(200);

		createBitimage(img, bitimg, val);

		// ImageJFunctions.show(bitimg).setTitle("KD-Tree input");

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());
		for (int i = 7; i < 8; ++i) {

			IterableInterval<BitType> view = Views.interval(bitimg, new long[] { i, i }, new long[] { 12 + i, 12 + i });
			list = getList(view);

			PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(n);

			PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(n);
			
			listonlyones = getvalueList(list, 1);
			listonlyzeros = getvalueList(list, 0);

		ArrayList<Point>	Xpointsort = getpointList(listonlyones);

		ArrayList<Point>	Ypointsort = getpointList(listonlyones);

			
			
			sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
			sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate

			System.out.println(Xpointsort);
			System.out.println(Ypointsort);
			final boolean biggeraxis = Xpointsort.size()-Ypointsort.size()>0;
			
			final int maxdir = biggeraxis? 0:1;
			
			Node<BitType> rootnode = makeNode(listonlyones, maxdir, 0);

			/*
			 * Cursor<BitType> cursor =
			 * rootnode.leftchild.leftchild.leftchild.leftchild.rightchild.
			 * Rightsublist.cursor(); while(cursor.hasNext()){ cursor.fwd();
			 * System.out.println(cursor.getDoublePosition(0));
			 * System.out.println(cursor.getDoublePosition(1));
			 * 
			 * }
			 */
			
			  Cursor<BitType> cursorleft = rootnode.Leftsublist.cursor();
			  System.out.println("MedianValue: "+rootnode.getMedianValue());
			  System.out.println("Direction: "+rootnode.getDirection());
			  
			  while(cursorleft.hasNext()){ cursorleft.fwd();
			  System.out.println("X-coordinate left:"
			  +cursorleft.getDoublePosition(0)); System.out.println(
			  "Y-coordinate left:" +cursorleft.getDoublePosition(1));
			  
			 } Cursor<BitType> cursorright = rootnode.Rightsublist.cursor();
			  
			  while(cursorright.hasNext()){ cursorright.fwd();
			  System.out.println("X-coordinate right:"
			  +cursorright.getDoublePosition(0)); System.out.println(
			  "Y-coordinate right:" +cursorright.getDoublePosition(1));
			  
			  }
			

			long startTimesec = System.currentTimeMillis();
			kdtree = TestConcisedistanceTransform(rootnode, list, imgout, new EucledianDistance());
			long endTimesec = System.currentTimeMillis();
			long totalTimesec = endTimesec - startTimesec;
			System.out.println(" O(nlog^2n) : " + totalTimesec);
			System.out.println(i);
			// new ImageJ();
			// ImageJFunctions.show(imgout).setTitle("KD-Tree output");

			long startTime = System.currentTimeMillis();
			bruteforce = TestBruteForce(list, listonlyones, brimgout, new EucledianDistance());
			long endTime = System.currentTimeMillis();
			long totalTime = endTime - startTime;
			System.out.println(" O(n^2) : " + totalTime);

			double count = 0;

			for (int index = 0; index < kdtree.size() - 1; ++index) {

				if (kdtree.get(index) - bruteforce.get(index) != 0) {
					count++;

				}
			}

			double rate = count / kdtree.size();

			System.out.println("Accuracy rate %: " + (1.0 - rate) * 100);
			/*
			 * PrintStream out = new PrintStream(new
			 * FileOutputStream("kdtree.txt")); System.setOut(out); for (int
			 * index =0; index<kdtree.size(); ++index)
			 * System.out.println(kdtree.get(index)); PrintStream outsec = new
			 * PrintStream(new FileOutputStream("brute.txt"));
			 * System.setOut(outsec); for (int indexsec =0;
			 * indexsec<kdtree.size(); ++indexsec)
			 * System.out.println(bruteforce.get(indexsec));
			 * 
			 */

			// ImageJFunctions.show(brimgout).setTitle("Brute-Tree output");

			/*
			 * 
			 * long startTimesec = System.currentTimeMillis(); kdtree =
			 * ConcisedistanceTransform(rootnode, listonlyzeros, imgout, new
			 * EucledianDistance()); long endTimesec =
			 * System.currentTimeMillis(); long totalTimesec = endTimesec -
			 * startTimesec; System.out.println(" O(nlog^2n) : " +
			 * totalTimesec); //new ImageJ();
			 * ImageJFunctions.show(imgout).setTitle("KD-Tree output");
			 * 
			 * long startTime = System.currentTimeMillis(); bruteforce =
			 * BruteForce(listonlyzeros, listonlyones, brimgout, new
			 * EucledianDistance()); long endTime = System.currentTimeMillis();
			 * long totalTime = endTime - startTime; System.out.println(
			 * " O(n^2) : " + totalTime);
			 * 
			 * ImageJFunctions.show(brimgout).setTitle("Brute-Tree output");
			 * 
			 * double count = 0;
			 * 
			 * for (int index = 0; index < kdtree.size() - 1; ++index) {
			 * 
			 * if (kdtree.get(index) - bruteforce.get(index) != 0) { count++;
			 * 
			 * } }
			 * 
			 * double rate = count / kdtree.size();
			 * 
			 * System.out.println("Accuracy rate %: " + (1.0 - rate) * 100);
			 */
		}
	}
}
