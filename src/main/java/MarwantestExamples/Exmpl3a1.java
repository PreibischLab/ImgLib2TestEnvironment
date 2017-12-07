package MarwantestExamples;

import java.util.Iterator;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.RealSum;

public class Exmpl3a1 {

	public <T extends RealType<T> & NativeType<T>> Exmpl3a1() throws ImgIOException {

		Img<T> img = (Img<T>) new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new ArrayImgFactory<FloatType>());
		T min = img.firstElement().createVariable();
		T max = img.firstElement().createVariable();

		computeMinMax(img,min,max);
		
		final double avg = computeAverage(img);
		
		System.out.println("Min:"+min);
		System.out.println("Max:"+max);
		System.out.println("Avg:"+avg);
		
		
	}
	private <T extends RealType<T>> double computeAverage(Img<T> input) {
		RealSum realSum = new RealSum();
		long count = 0;
		for (final T type:input) {
			realSum.add(type.getRealDouble());
			++count;
		}
		return realSum.getSum()/count;
	}
	public <T extends Comparable<T> & Type<T>> void computeMinMax(Img<T> input, T min, T max) {

		
		final Iterator<T> iterator = input.iterator();
		
		T type = iterator.next();
		
		min.set(type);
		max.set(type);
		
		while(iterator.hasNext()) {
			type= iterator.next();
			if (type.compareTo(min)<0) {
				min.set(type);
			}
			if (type.compareTo(max)>0) {
				max.set(type);
			}
		}
	}
	public static void main(String[] args) throws ImgIOException {
		new ImageJ();
		new Exmpl3a1();
	}
}
