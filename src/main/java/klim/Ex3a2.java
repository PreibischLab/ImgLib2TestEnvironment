package klim;

import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;

public class Ex3a2 {

	public Ex3a2(){
		ArrayList<FloatType> list = new ArrayList<FloatType>();
		for (int i = 0; i < 10; ++i)
			list.add(new FloatType(i));
		
		FloatType minVal = new FloatType();
		FloatType maxVal = new FloatType();
		computeMinMax(list, minVal, maxVal);
		System.out.println(minVal + " <- min");
		System.out.println(maxVal + " <- max");
	}
	
	public <T extends Comparable<T> & Type<T>> void computeMinMax(final Iterable<T> input, final T min, final T max){
		final Iterator<T> iterator = input.iterator();
		T type = iterator.next();
		min.set(type);
		max.set(type);
		
		while(iterator.hasNext()){
			type = iterator.next();
			if (type.compareTo(max) > 0){
				max.set(type);
			}
			if (type.compareTo(min) < 0){
				min.set(type);
			}	
		}
	
		
	} 
	
	public static void main(String[] args){
		new Ex3a2();
	}
}
