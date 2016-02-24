package varun;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import varun.MyKDtree.Distance;
import varun.MyKDtree.Node;

public class KDtreejunk {
	
	/*

	/******
	 * Returns a root tree, I do this to initialize an ArrayList<Node<T>> in the
	 * main program which I overwrite later to include all the subtrees
	 ******/
/*
	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction) {

		int n = list.numDimensions();

		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
	/*
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

			Cursor<T> rightTreecursor = RightTree.localizingCursor();
			rightTreecursor.fwd();

			Node<T> node = new Node<T>(pivotElement, rightTreecursor, direction, LeftTree, RightTree);

			return node;

		}

	}

	

	/**********
	 * Starting the distance transform routine, done in a sort of non-Java
	 * fashion, by getting a needed object out of a void method by creating an
	 * empty ArrayList of the object and then overwritting it in the program.
	 * For a better tasting implementation see the overloaded version below
	 * which creates a NN object and uses getters and setters to get the correct
	 * distance transform.
	 * 
	 * @throws FileNotFoundException
	 **********/
/*
	public static <T extends RealType<T>> void distanceTransform(PointSampleList<BitType> list,
			PointSampleList<BitType> listzerosorones, RandomAccessibleInterval<T> imgout, final Distance dist)
					throws FileNotFoundException {

		Node<BitType> rootnode;

		rootnode = makeNode(list, 0);
		PrintStream out = new PrintStream(new FileOutputStream("KDtreemindist.txt"));
		System.setOut(out);

		double Bestdistsquared = Double.MAX_VALUE;

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		double distance = 0;
		while (zerooronelistcursor.hasNext()) {
			zerooronelistcursor.fwd();
			double mindistance = Double.MAX_VALUE;
			double farmindistance = Double.MAX_VALUE;
			outbound.setPosition(zerooronelistcursor);

			ArrayList<Node<BitType>> nearnodelist = new ArrayList<Node<BitType>>(4);
			ArrayList<Node<BitType>> farnodelist = new ArrayList<Node<BitType>>(4);

			closestNode(zerooronelistcursor, rootnode, nearnodelist, Bestdistsquared);
			furtherNode(zerooronelistcursor, rootnode, farnodelist, Bestdistsquared);
			PointSampleList<BitType> singletree = combineTrees(nearnodelist.get(0));
			PointSampleList<BitType> singlefartree;
			if (farnodelist != null)
				singlefartree = combineTrees(farnodelist.get(0));
			else
				singlefartree = null;

			Cursor<BitType> singlecursor = singletree.cursor();

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerooronelistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			if (singlefartree != null) {

				Cursor<BitType> singlefarcursor = singlefartree.cursor();
				double fardistance = 0;
				while (singlefarcursor.hasNext()) {
					singlefarcursor.fwd();

					fardistance = dist.getDistance(zerooronelistcursor, singlefarcursor);
					farmindistance = Math.min(fardistance, farmindistance);

				}

			}
			final double actualmindistance = Math.min(mindistance, farmindistance);

			outbound.get().setReal(actualmindistance);

		}

		final Cursor<BitType> listcursor = list.localizingCursor();
		while (listcursor.hasNext()) {
			listcursor.fwd();
			outbound.setPosition(listcursor);
			outbound.get().setReal(0);

		}

	}
	/*
	/***********
	 * Returns the node closest to the testpoint stores as a single entry in the
	 * arrayList, done in a sort of "non-java" fashion
	 ***********/
/*
	public static <T extends RealType<T>> void closestNode(RealLocalizable testpoint, Node<T> Trees,
			ArrayList<Node<T>> list, double Bestdistsquared) {

		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else
			otherdirection = direction + 1;

		double locationdiff = (testpoint.getDoublePosition(direction) - Trees.nodePoint[direction]);

		double dist = 0;
		for (int d = 0; d < n; ++d) {

			dist += (testpoint.getDoublePosition(d) - Trees.nodePoint[d])
					* (testpoint.getDoublePosition(d) - Trees.nodePoint[d]);
		}

		final Node<T> finalnode;

		final double mindistsquared = Math.min(dist, Bestdistsquared);

		if (dist < Bestdistsquared) {

			finalnode = Trees;
			Bestdistsquared = mindistsquared;

		}

		else
			finalnode = null;

		final boolean leftbranchsearch = locationdiff < 0;

		final PointSampleList<T> searchBranch = leftbranchsearch ? Trees.LeftTree : Trees.RightTree;

		Node<T> nearnode;
		if ((searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1) > 2) {

			nearnode = makeNode(searchBranch, otherdirection);
			if (nearnode != null)
				closestNode(testpoint, nearnode, list, mindistsquared);

		}

		if (finalnode != null)
			list.add(finalnode);

	}

	public static <T extends RealType<T>> void furtherNode(RealLocalizable testpoint, Node<T> Trees,
			ArrayList<Node<T>> list, double Bestdistsquared) {
		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else
			otherdirection = direction + 1;

		double locationdiff = (testpoint.getDoublePosition(direction) - Trees.nodePoint[direction]);

		double axisdiff = locationdiff * locationdiff;

		double dist = 0;
		for (int d = 0; d < n; ++d) {

			dist += (testpoint.getDoublePosition(d) - Trees.nodePoint[d])
					* (testpoint.getDoublePosition(d) - Trees.nodePoint[d]);
		}

		final Node<T> finalnode;

		final double mindistsquared = Math.min(dist, Bestdistsquared);

		if (dist < Bestdistsquared) {

			finalnode = Trees;

		}

		else
			finalnode = null;

		final boolean leftbranchsearch = locationdiff < 0;

		final PointSampleList<T> nonsearchBranch = leftbranchsearch ? Trees.RightTree : Trees.LeftTree;

		Node<T> farnode;

		if (axisdiff <= mindistsquared
				&& (nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1) > 2) {
			farnode = makeNode(nonsearchBranch, otherdirection);
			if (farnode != null)
				furtherNode(testpoint, farnode, list, mindistsquared);
		}

		if (finalnode != null)
			list.add(finalnode);

	}
*/
}
