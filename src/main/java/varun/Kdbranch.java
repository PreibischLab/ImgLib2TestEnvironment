package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import net.imglib2.PointSampleList;
import net.imglib2.RealPointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Kdbranch {
	
	
	private Kdnode defaultNode = null;
	
	private int k=0;
	
	
	
	
	
	public static class Kdnode{
		
		private final Point p;
		
		private final int direction;
		
		public  Kdnode(Point p){
			this.p=p;
			this.direction=0;
			}
		
		public  Kdnode(Point p, int direction){
			this.p=p;
			this.direction=direction;
			}
		
		
		
		
	}
	
	
	// Default Constructor
	public Kdbranch(){
		
	}
	
	// Constructor for creating tree from a PointSampleList.
	 public  <T extends RealType<T>>   Kdbranch(PointSampleList<T> list) {
	        super();
	        defaultNode = createNode(list,k);
	    }
	 
	// To Create node from a list of points IterableInterval<T> or I would probably change this to PointSampleList<T>
	public static <T extends RealType<T>> Kdnode createNode(IterableInterval<T> list, int direction) {
		
		long listSize=list.size();
		
		
PointSampleList<T> childA= new PointSampleList<T>((int) listSize);

PointSampleList<T> childB= new PointSampleList<T>((int) listSize);

		Kdnode node=null;
		
		

		int meanIndex = (int) list.size() / 2;
		
		Point meanList= new Point(meanIndex);
		
		node=new Kdnode(meanList, direction); // Create a node at the mean value
		

		 if ((meanIndex - 1) >= 0 && childA.size() > 0) {

		 createNode(childA,direction+1);

	}

	 if ((meanIndex + 1) <= (list.size() - 1) && childB.size() > 0) {
	 createNode(childB,direction + 1);
	 }

	 
	 return node;
	 
	 }

	// How to create a pixel value list from an IterableInterval<T>?

	public static <T extends RealType<T>> void getSortedImage(RandomAccessibleInterval<T> img,
			RandomAccessibleInterval<T> imgout) {

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
		
		
		
		final RandomAccess<T> outbound = imgout.randomAccess();

		first.fwd();
		
		
		
		
		
		outbound.setPosition(first);

		while (first.hasNext()) {
			first.fwd();
			
		
			
			
			
				
			
			if (outbound.get().compareTo(first.get()) >= 0)
				
				outbound.setPosition(first);
		

			outbound.get().set(first.get());
			
			
		
		}
		
	}

	public static void main(String[] args) {
		
		
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bud.jpg"));

		Img<FloatType> imgout = new CellImgFactory<FloatType>().create(img, new FloatType());

		ImageJFunctions.show(img).setTitle("Original_Image");

		getSortedImage(img, imgout);

		ImageJFunctions.show(imgout).setTitle("Permuted_Image");

	}

}