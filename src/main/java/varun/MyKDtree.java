package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import com.sun.tools.javac.api.Formattable.LocalizedString;

import ij.ImageJ;
import net.imglib2.PointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTreeNode;
import net.imglib2.Localizable;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPointSampleList;
import net.imglib2.Sampler;
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

	

	public static <T extends RealType<T>> void split(PointSampleList<T> list, int direction) {

		int n = list.numDimensions();

		if (list.size() <= 1)
			return;

		else {

			/****
			 * To ward against running over the dimensionality, creating some
			 * local restrictions on the global variable direction
			 ****/
			if (direction == list.numDimensions())
				direction = 0;

			// the first element belonging to the right list childB
			final int splitIndex = (int) list.size() / 2;

			final PointSampleList<T> childA = new PointSampleList<T>(n);
			final PointSampleList<T> childB = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			int index = 0;
			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

				if (index < splitIndex) {

					childA.add(cord, listCursor.get().copy());

					// System.out.println("childA: "+listCursor.get());

				} else

				{

					childB.add(cord, listCursor.get().copy());
					// System.out.println("childB: "+listCursor.get());
				}
				index++;
			}

			split(childA, direction);

			split(childB, direction);

			mergeList(list, childA, childB, direction);
		}

	}

	/// ***** Returns a sorted list *********////
	public static <T extends RealType<T>> void mergeList(PointSampleList<T> list, PointSampleList<T> listA,
			PointSampleList<T> listB, int direction) {

		final Cursor<T> cursorA = listA.localizingCursor();
		final Cursor<T> cursorB = listB.localizingCursor();
		final Cursor<T> cursor = list.localizingCursor();

		cursorA.fwd();
		cursorB.fwd();

		boolean cannotMoveOn = false;

		do {
			// here is where you decide what you sort after
			if (cursorA.getDoublePosition(direction) < (cursorB.getDoublePosition(direction))) {

				cursor.fwd();
				cursor.get().set(cursorA.get());
				if (cursorA.hasNext())
					cursorA.fwd();
				else {
					cannotMoveOn = true;

					// move cursorB until the end
					boolean stopped = false;
					do {
						cursor.fwd();
						cursor.get().set(cursorB.get());
						if (cursorB.hasNext())
							cursorB.fwd();
						else
							stopped = true;
					} while (stopped == false);
				}

			}

			else

			{

				cursor.fwd();
				cursor.get().set(cursorB.get());
				if (cursorB.hasNext())
					cursorB.fwd();
				else {
					cannotMoveOn = true;

					// move cursorA until the end
					boolean stopped = false;
					do {
						cursor.fwd();
						cursor.get().set(cursorA.get());
						if (cursorA.hasNext())
							cursorA.fwd();
						else
							stopped = true;
					} while (stopped == false);
				}

			}

		} while (cannotMoveOn == false);
	}

	

	public static class Node<T> {

		public final int n;

		public final double medianValue;

		public final double[] nodePoint;

		public final int direction;

		public final PointSampleList<T> LeftTree;

		public final PointSampleList<T> RightTree;

		public Node(final double medianValue, RealLocalizable point, final int direction,
				final PointSampleList<T> LeftTree, final PointSampleList<T> RightTree) {
			assert LeftTree.numDimensions() == RightTree.numDimensions();
			this.n = LeftTree.numDimensions();
			this.nodePoint = new double[n];
			point.localize(nodePoint);

			this.medianValue = medianValue;
			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;
		}

		public Node() {
			super();
			this.n = 0;
			this.medianValue = 0;
			this.nodePoint = new double[n];
			this.direction = 0;
			this.LeftTree = null;
			this.RightTree = null;
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

		public void localize(double[] point) {

			for (int d = 0; d < n; ++d)
				point[d] = (double) nodePoint[d];
		}

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	public static class searchNode<T> {

		private final int n;

		public final double[] Position;

		protected PointSampleList<T> list;

		protected Node<T> finalnode;

		protected double Bestdistsquared;

		protected Node<T> farfinalnode;

		protected double Bestfardistsquared;

		public searchNode(final PointSampleList<T> list) {
			n = list.numDimensions();
			Position = new double[n];
			this.list = list;

		}

		public int getNumdimensions() {
			return n;
		}

		public void search(final RealLocalizable cursor) {
			cursor.localize(Position);
			Bestdistsquared = Double.MAX_VALUE;
			closestNode(makeNode(list, 0));
		}

		public void searchfar(final RealLocalizable cursor) {
			cursor.localize(Position);
			Bestfardistsquared = Double.MAX_VALUE;
			furtherNode(makeNode(list, 0));
		}

		public double getBestdist() {
			return Bestdistsquared;
		}

		public Node<T> getfinalnode() {
			return finalnode;
		}

		public double getfarBestdist() {
			return Bestfardistsquared;
		}

		public Node<T> getfarfinalnode() {
			return farfinalnode;
		}

		private void closestNode(final Node<T> currentBest) {
			int direction = currentBest.direction;

			int n = currentBest.getnumDimensions();
			int otherdirection;

			if (direction == n - 1)

				otherdirection = 0;

			else
				otherdirection = direction + 1;

			double dist = 0;
			for (int d = 0; d < n; ++d) {

				dist += (Position[d] - currentBest.nodePoint[d]) * (Position[d] - currentBest.nodePoint[d]);
			}

			if (dist < Bestdistsquared) {
				Bestdistsquared = dist;
				finalnode = currentBest;
			}

			final double locationdiff = Position[currentBest.direction] - currentBest.nodePoint[currentBest.direction];
			final double axisdiff = locationdiff * locationdiff;
			final boolean leftbranchsearch = locationdiff < 0;

			
			final PointSampleList<T> searchBranch = leftbranchsearch ? currentBest.LeftTree : currentBest.RightTree;
			final PointSampleList<T> nonsearchBranch = leftbranchsearch ? currentBest.RightTree : currentBest.LeftTree;
			Node<T> nearnode;
			if ((searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1) > 2) {
				nearnode = makeNode(searchBranch, otherdirection);
				if (nearnode != null)
					closestNode(nearnode);

			}
			

		}

		private void furtherNode(final Node<T> currentBest) {
			int direction = currentBest.direction;

			int n = currentBest.getnumDimensions();
			int otherdirection;

			if (direction == n - 1)

				otherdirection = 0;

			else
				otherdirection = direction + 1;

			double dist = 0;
			for (int d = 0; d < n; ++d) {

				dist += (Position[d] - currentBest.nodePoint[d]) * (Position[d] - currentBest.nodePoint[d]);
			}

			if (dist < Bestfardistsquared) {
				Bestfardistsquared = dist;
				farfinalnode = currentBest;
			}

			final double locationdiff = Position[currentBest.direction] - currentBest.nodePoint[currentBest.direction];
			final double axisdiff = locationdiff * locationdiff;
			final boolean leftbranchsearch = locationdiff < 0;

			// search the far branch

			final PointSampleList<T> nonsearchBranch = leftbranchsearch ? currentBest.RightTree : currentBest.LeftTree;

			Node<T> farnode;
			if (axisdiff <= Bestfardistsquared
					&& (nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1) > 2) {
				farnode = makeNode(nonsearchBranch, otherdirection);
				if (farnode != null)
					closestNode(farnode);
			}
		}

		private double getMedian(PointSampleList<T> list, int direction) {

			final Cursor<T> listCursor = list.localizingCursor();

			final ArrayList<Double> values = new ArrayList<Double>();

			while (listCursor.hasNext()) {
				listCursor.fwd();

				values.add(listCursor.getDoublePosition(direction));

			}

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

		private Node<T> makeNode(PointSampleList<T> list, int direction) {

			int n = list.numDimensions();

			/****
			 * To ward against running over the dimensionality, creating some
			 * local restrictions on the global variable direction
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

						LeftTree.add(cord, listCursor.get());

					else

						RightTree.add(cord, listCursor.get());

				}

				Cursor<T> rightTreecursor = RightTree.localizingCursor();
				rightTreecursor.fwd();

				Node<T> node = new Node<T>(pivotElement, rightTreecursor, direction, LeftTree, RightTree);

				return node;

			}

		}

	}

	
	

	/********
	 * For a Node<T>, returns a single PointSampleList by combining the Left and
	 * Right Tree pairs into one
	 ***********/

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

			singleTree.add(treepointB, treecursorB.get().copy());

		}

		return singleTree;

	}


	/****
	 * This is more of a "Java" way of getting NN by using an NN object and
	 * using getters and setters to get the Best Node point
	 ****/
	public static <T extends RealType<T>> void ConcisedistanceTransform(PointSampleList<BitType> list,
			PointSampleList<BitType> listzerosorones, RandomAccessibleInterval<T> imgout, final Distance dist)
					throws FileNotFoundException {

	//	PrintStream out = new PrintStream(new FileOutputStream("conKDtreemindist.txt"));
	//	System.setOut(out);

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		final searchNode<BitType> Bestnode = new searchNode<BitType>(list);

		final searchNode<BitType> farBestnode = new searchNode<BitType>(list);

		double distance = 0;
		double fardistance = 0;
		while (zerooronelistcursor.hasNext()) {
			zerooronelistcursor.fwd();
			double mindistance = Double.MAX_VALUE;
			double farmindistance = Double.MAX_VALUE;

			outbound.setPosition(zerooronelistcursor);

			Bestnode.search(zerooronelistcursor);

			PointSampleList<BitType> singletree = combineTrees(Bestnode.getfinalnode());

		
			
			Cursor<BitType> singlecursor = singletree.cursor();

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerooronelistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			farBestnode.searchfar(zerooronelistcursor);

			PointSampleList<BitType> farsingletree = combineTrees(farBestnode.getfarfinalnode());

			Cursor<BitType> farsinglecursor = farsingletree.cursor();

			while (farsinglecursor.hasNext()) {
				farsinglecursor.fwd();

				fardistance = dist.getDistance(zerooronelistcursor, farsinglecursor);
				farmindistance = Math.min(fardistance, farmindistance);

			}

			final double actualmindistance = Math.min(farmindistance, mindistance);

			

			outbound.get().setReal(actualmindistance);

		}

		final Cursor<BitType> listcursor = list.localizingCursor();
		while (listcursor.hasNext()) {
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
		
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());

		FloatType val = new FloatType(200);

		ImageJFunctions.show(img).setTitle("KD-Tree input");
		
		createBitimage(img, bitimg, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		IterableInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 100, 100 });

		list = getList(view);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(bitimg.numDimensions());

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(bitimg.numDimensions());

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		split(listonlyones, 0);
		split(listonlyones, 1);
		 long startTime = System.currentTimeMillis();
		ConcisedistanceTransform(listonlyones, listonlyzeros, imgout, new EucledianDistance());

	//	new ImageJ();

	

		ImageJFunctions.show(imgout).setTitle("KD-Tree output");
		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println(totalTime);
		
		
	}
}