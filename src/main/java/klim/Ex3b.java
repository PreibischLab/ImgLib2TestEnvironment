package klim;

import java.io.File;

import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.RealSum;
import util.ImgLib2Util;

public class Ex3b {

	public <T extends RealType<T> & NativeType <T>> Ex3b(){
		Img<T> img = (Img<T>) ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final double avg = computeAverage(img);
		
		System.out.println("average " + avg);
	}
	
	public <T extends RealType<T>> double computeAverage(final Iterable<T> input){
		final RealSum realSum = new RealSum();
		long count  = 0;
		
		for (final T type : input){
			count++;
			realSum.add(type.getRealDouble());
		}
		
		return realSum.getSum()/count;
	}
	
	public static void main(String[] args){
		new Ex3b();
	}
}
