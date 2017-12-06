package marwan;

import java.awt.Rectangle;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.real.FloatType;

public class Portion extends Object {

	private RandomAccessibleInterval<FloatType> view;
	private Rectangle shape;
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
	
	
	
}

