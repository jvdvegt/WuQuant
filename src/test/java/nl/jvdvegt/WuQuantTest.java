package nl.jvdvegt;

import static java.awt.image.BufferedImage.TYPE_INT_RGB;

import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.apache.commons.imaging.formats.png.PngImagingParameters;
import org.apache.commons.imaging.formats.png.PngWriter;
import org.apache.commons.imaging.palette.Palette;
import org.apache.commons.imaging.palette.PaletteFactory;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class WuQuantTest {
    @Test
    void test() throws IOException {
        // Note: quantization is extremely hard for this image,
        // and it does not show in no way that Wu's algorithm is better.
        //
        // This test is merely here to show how to use WuQuant.java.
        final BufferedImage image = createDummyImage();

        final byte[] normalBytes = createPngBytes(image, new PaletteFactory());
        final byte[] wuBytes = createPngBytes(image, new WuPaletteFactory());

//        Files.write(Paths.get("normal.png"), normalBytes);
//        Files.write(Paths.get("wu.png"), wuBytes);
        Assertions.assertTrue(normalBytes.length >  wuBytes.length);
    }

    static final class WuPaletteFactory extends PaletteFactory {
        @Override
        public Palette makeQuantizedRgbPalette(final BufferedImage src, final int max) {
            return new WuQuant(src).createPalette(max);
        }
    }

    private static byte[] createPngBytes(final BufferedImage image, final PaletteFactory factory) throws IOException {
        try (ByteArrayOutputStream os = new ByteArrayOutputStream()) {
            final PngImagingParameters params = new PngImagingParameters();
            params.setForceIndexedColor(true);
            new PngWriter().writeImage(image, os, params, factory);

            return os.toByteArray();
        }
    }

    private static BufferedImage createDummyImage() {
        final BufferedImage image = new BufferedImage(256, 256, TYPE_INT_RGB);
        for (int i = 0; i < image.getWidth(); i++) {
            for (int j = 0; j < image.getHeight(); j++) {
                float h = i / (float) image.getWidth();
                float s = (image.getHeight() - j) / (float) image.getHeight();
                float l = 0.5f;
                int rgb = hslToRgb(h, s, l);
                image.setRGB(i, j, rgb);
            }
        }
        return image;
    }

    // From https://stackoverflow.com/a/29316972 :
    private static int hslToRgb(float h, float s, float l){
        float r, g, b;

        if (s == 0f) {
            r = g = b = l; // achromatic
        } else {
            float q = l < 0.5f ? l * (1 + s) : l + s - l * s;
            float p = 2 * l - q;
            r = hueToRgb(p, q, h + 1f/3f);
            g = hueToRgb(p, q, h);
            b = hueToRgb(p, q, h - 1f/3f);
        }
        int rgb = to255(r) << 16 | to255(g) << 8 | to255(b);
        return rgb;
    }

    private static int to255(float v) {
        return (int) Math.min(255, 256 * v);
    }

    /** Helper method that converts hue to rgb */
    private static float hueToRgb(float p, float q, float t) {
        if (t < 0f)
            t += 1f;
        if (t > 1f)
            t -= 1f;
        if (t < 1f/6f)
            return p + (q - p) * 6f * t;
        if (t < 1f/2f)
            return q;
        if (t < 2f/3f)
            return p + (q - p) * (2f/3f - t) * 6f;
        return p;
    }
}
