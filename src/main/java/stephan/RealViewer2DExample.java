/*
 * #%L
 * ImgLib2: a general-purpose, multidimensional image processing library.
 * %%
 * Copyright (C) 2009 - 2015 Tobias Pietzsch, Stephan Preibisch, Barry DeZonia,
 * Stephan Saalfeld, Curtis Rueden, Albert Cardona, Christian Dietz, Jean-Yves
 * Tinevez, Johannes Schindelin, Jonathan Hale, Lee Kamentsky, Larry Lindsey, Mark
 * Hiner, Michael Zinsmaier, Martin Horn, Grant Harris, Aivar Grislis, John
 * Bogovic, Steffen Jaensch, Stefan Helfrich, Jan Funke, Nick Perry, Mark Longair,
 * Melissa Linkert and Dimiter Prodanov.
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

package stephan;
import java.io.File;

import net.imglib2.RandomAccessible;
import net.imglib2.RealRandomAccessible;
import net.imglib2.converter.RealARGBConverter;
import net.imglib2.img.Img;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.ui.overlay.LogoPainter;
import net.imglib2.ui.viewer.InteractiveRealViewer2D;
import net.imglib2.view.Views;
import util.ImgLib2Util;

/**
 * Adapted from MandelbrotRealViewer2DExample in imglib2-tutorials
 * 
 * @author spreibi
 *
 * @param <T>
 */
public class RealViewer2DExample< T extends NumericType< T > & NativeType< T > >
{
	final static public void main( final String[] args )
	{
		final int width = 800;
		final int height = 600;

		final Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/bridge.png" ) );
		final RandomAccessible< FloatType > ra = Views.extendMirrorSingle( img );
		final RealRandomAccessible< FloatType > rra = Views.interpolate( ra, new NLinearInterpolatorFactory< FloatType >() );

		final AffineTransform2D transform = new AffineTransform2D();
		transform.scale( 200 );
		transform.translate( width / 2.0, height / 2.0 );

		final RealARGBConverter< FloatType > converter = new RealARGBConverter< FloatType >( 0, 255 );
		final InteractiveRealViewer2D< FloatType > viewer = new InteractiveRealViewer2D< FloatType >( width, height, rra, transform, converter );
		viewer.getDisplayCanvas().addOverlayRenderer( new LogoPainter() );
		viewer.requestRepaint();
	}
}
