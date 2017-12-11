package marwan;

import java.awt.Rectangle;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;

public class Portion extends Object {

	private RandomAccessibleInterval<FloatType> view;
	private Rectangle shape;
	private int slice;
	private int dimenssion;
	
	

	public Portion(RandomAccessibleInterval<FloatType> view, Rectangle shape,int dimenssion, int slice) {
		super();
		this.view = view;
		this.shape = shape;
		this.slice = slice;
		this.dimenssion = dimenssion;
	}
	public Portion(RandomAccessibleInterval<FloatType> view, Rectangle shape) {
		super();
		this.view = view;
		this.shape = shape;
	}
	public RandomAccessibleInterval<FloatType> getView() {
		return view;
	}
	public void setView(RandomAccessibleInterval<FloatType> view) {
		this.view = view;
	}
	public Rectangle getShape() {
		return shape;
	}
	public void setShape(Rectangle shape) {
		this.shape = shape;
	}
	
	public int getSlice() {
		return slice;
	}
	public void setSlice(int slice) {
		this.slice = slice;
	}
	public int getDimenssion() {
		return dimenssion;
	}
	public void setDimenssion(int dimenssion) {
		this.dimenssion = dimenssion;
	}
	
}

