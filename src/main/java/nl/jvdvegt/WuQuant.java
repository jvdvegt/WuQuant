package nl.jvdvegt;

import java.awt.image.BufferedImage;

import org.apache.commons.imaging.palette.Palette;

/**
 * Java version of Xiaolin Wu's color Quantization algorithm. Code is cleaned a bit but resembles the
 * <a href="https://gist.github.com/bert/1192520">original C code</a>.
 * <p>
 * Free to distribute, comments and suggestions are appreciated.
 * </p>
 * Original header follows.
 */

/*
        C Implementation of Wu's Color Quantizer (v. 2)
        (see Graphics Gems vol. II, pp. 126-133)

Author:    Xiaolin Wu
           Dept. of Computer Science
           Univ. of Western Ontario
           London, Ontario N6A 5B7
           wu@csd.uwo.ca

Algorithm: Greedy orthogonal bipartition of RGB space for variance
           minimization aided by inclusion-exclusion tricks.
           For speed no nearest neighbor search is done. Slightly
           better performance can be expected by more sophisticated
           but more expensive versions.

The author thanks Tom Lane at Tom_Lane@G.GP.CS.CMU.EDU for much of
additional documentation and a cure to a previous bug.
Free to distribute, comments and suggestions are appreciated.
*/
public final class WuQuant {
    private static final int RED = 2;
    private static final int GREEN = 1;
    private static final int BLUE = 0;

    private static final int SIZE = 33;

    private final int transparentColor; // Transparent color

    /* Histogram is in elements 1..HISTSIZE along each axis,
     * element 0 is for base or marginal value
     * NB: these must start out 0!
     */
    private final float[][][] m2 = new float[SIZE][SIZE][SIZE];
    private final long[][][] wt = new long[SIZE][SIZE][SIZE];
    private final long[][][] mr = new long[SIZE][SIZE][SIZE];
    private final long[][][] mg = new long[SIZE][SIZE][SIZE];
    private final long[][][] mb = new long[SIZE][SIZE][SIZE];

    public WuQuant(final BufferedImage src) {
        // First we build a histogram
        transparentColor = hist3D(src);

        /* We now convert histogram into moments so that we can rapidly calculate
         * the sums of the above quantities over any desired box.
         */
        m3D();
    }

    public Palette createPalette(final int paletteSize) {
        final Box[] cube = new Box[paletteSize];
        for (int i = 0; i < paletteSize; i++) {
            cube[i] = new Box();
        }

        cube[0].r0 = cube[0].g0 = cube[0].b0 = 0;
        cube[0].r1 = cube[0].g1 = cube[0].b1 = SIZE - 1;

        final int resultSize = transparentColor >= 0 ? paletteSize - 1 : paletteSize;
        final float[] vv = new float[paletteSize];
        int next = 0;
        for (int i = 1; i < resultSize; ++i) {
            if (cut(cube[next], cube[i])) {
                /* volume test ensures we won't try to cut one-cell box */
                vv[next] = (cube[next].vol > 1) ? variance(cube[next]) : 0.0f;
                vv[i] = (cube[i].vol > 1) ? variance(cube[i]) : 0.0f;
            }
            else {
                vv[next] = 0.0f;   /* don't try to split this box again */
                i--;              /* didn't create box i */
            }
            next = 0;
            float temp = vv[0];
            for (int k = 1; k <= i; ++k) {
                if (vv[k] > temp) {
                    temp = vv[k];
                    next = k;
                }
            }
            if (temp <= 0.0f) {
                break;
            }
        }

        final int[] tag = new int[SIZE * SIZE * SIZE];

        final int[] lut = new int[resultSize];
        for (int k = 0; k < resultSize; ++k) {
            mark(cube[k], k, tag);

            final long weight = volume(cube[k], wt);
            if (weight > 0) {
                final int lutR = (int) (volume(cube[k], mr) / weight);
                final int lutG = (int) (volume(cube[k], mg) / weight);
                final int lutB = (int) (volume(cube[k], mb) / weight);
                lut[k] = (255 << 24) | (lutR << 16) | (lutG << 8) | lutB;
            }
            else {
                lut[k] = 0;
            }
        }

        if (transparentColor >= 0) {
            lut[paletteSize] = transparentColor; // Set the transparent color
        }
        return new Palette() {
            @Override
            public int getPaletteIndex(final int rgb) {
                final int r = ((rgb >> 16) & 0xff);
                final int g = ((rgb >> 8) & 0xff);
                final int b = (rgb & 0xff);
                final int inr = (r >> 3) + 1;
                final int ing = (g >> 3) + 1;
                final int inb = (b >> 3) + 1;
                final int tagIndex = (inr << 10) + (inr << 6) + inr + (ing << 5) + ing + inb;
                return tag[tagIndex];
            }

            @Override
            public int getEntry(final int index) {
                return lut[index];
            }

            @Override
            public int length() {
                return lut.length;
            }
        };
    }

    /* Histogram is in elements 1..HISTSIZE along each axis,
     * element 0 is for base or marginal value
     * NB: these must start out 0!
     */
    private int hist3D(final BufferedImage image) {
        /* build 3-D color histogram of counts, r/g/b, c^2 */
        final int[] table = new int[256];

        for (int i = 0; i < 256; ++i) {
            table[i] = i * i;
        }
        int transparentColor = -1;
        final int width = image.getWidth();
        final int height = image.getHeight();

        final int[] row = new int[width];

        for (int y = 0; y < height; y++) {
            image.getRGB(0, y, width, 1, row, 0, width);
            for (int x = 0; x < width; x++) {
                final int argb = row[x];
                if ((argb >>> 24) < 0x80) { // Transparent
                    if (transparentColor < 0) {   // Find the transparent color
                        transparentColor = argb;
                    }
                }
                final int r = ((argb >> 16) & 0xff);
                final int g = ((argb >> 8) & 0xff);
                final int b = (argb & 0xff);
                final int inr = (r >> 3) + 1;
                final int ing = (g >> 3) + 1;
                final int inb = (b >> 3) + 1;

                /*[inr][ing][inb]*/
                ++wt[inr][ing][inb];
                mr[inr][ing][inb] += r;
                mg[inr][ing][inb] += g;
                mb[inr][ing][inb] += b;
                m2[inr][ing][inb] += table[r] + table[g] + table[b];
            }
        }
        return transparentColor;
    }

    /* At conclusion of the histogram step, we can interpret
     *   wt[r][g][b] = sum over voxel of P(c)
     *   mr[r][g][b] = sum over voxel of r*P(c)  ,  similarly for mg, mb
     *   m2[r][g][b] = sum over voxel of c^2*P(c)
     * Actually each of these should be divided by 'size' to give the usual
     * interpretation of P() as ranging from 0 to 1, but we needn't do that here.
     */
    private void m3D() {
        /* compute cumulative moments. */
        for (int r = 1; r < SIZE; ++r) {
            final int[] area = new int[SIZE];
            final int[] areaR = new int[SIZE];
            final int[] areaG = new int[SIZE];
            final int[] areaB = new int[SIZE];
            final float[] area2 = new float[SIZE];

            for (int g = 1; g < SIZE; ++g) {
                float line2 = 0;
                long line = 0;
                long lineR = 0;
                long lineG = 0;
                long lineB = 0;
                for (int b = 1; b < SIZE; ++b) {
                    line += wt[r][g][b];
                    lineR += mr[r][g][b];
                    lineG += mg[r][g][b];
                    lineB += mb[r][g][b];
                    line2 += m2[r][g][b];

                    area[b] += line;
                    areaR[b] += lineR;
                    areaG[b] += lineG;
                    areaB[b] += lineB;
                    area2[b] += line2;

                    wt[r][g][b] = wt[r - 1][g][b] + area[b];
                    mr[r][g][b] = mr[r - 1][g][b] + areaR[b];
                    mg[r][g][b] = mg[r - 1][g][b] + areaG[b];
                    mb[r][g][b] = mb[r - 1][g][b] + areaB[b];
                    m2[r][g][b] = m2[r - 1][g][b] + area2[b];
                }
            }
        }
    }

    private static long volume(final Box cube, final long[][][] mmt) {
        /* Compute sum over a box of any given statistic */
        return (mmt[cube.r1][cube.g1][cube.b1]
            - mmt[cube.r1][cube.g1][cube.b0]
            - mmt[cube.r1][cube.g0][cube.b1]
            + mmt[cube.r1][cube.g0][cube.b0]
            - mmt[cube.r0][cube.g1][cube.b1]
            + mmt[cube.r0][cube.g1][cube.b0]
            + mmt[cube.r0][cube.g0][cube.b1]
            - mmt[cube.r0][cube.g0][cube.b0]);
    }

    /* The next two routines allow a slightly more efficient calculation
     * of Vol() for a proposed subbox of a given box.  The sum of Top()
     * and Bottom() is the Vol() of a subbox split in the given direction
     * and with the specified new upper bound.
     */

    private static long bottom(final Box cube, final int color, final long[][][] mmt) {
        /* Compute part of Vol(cube, mmt) that doesn't depend on r1, g1, or b1 */
        /* (depending on color) */
        switch (color) {
            case RED:
                return (-mmt[cube.r0][cube.g1][cube.b1]
                    + mmt[cube.r0][cube.g1][cube.b0]
                    + mmt[cube.r0][cube.g0][cube.b1]
                    - mmt[cube.r0][cube.g0][cube.b0]);
            case GREEN:
                return (-mmt[cube.r1][cube.g0][cube.b1]
                    + mmt[cube.r1][cube.g0][cube.b0]
                    + mmt[cube.r0][cube.g0][cube.b1]
                    - mmt[cube.r0][cube.g0][cube.b0]);
            case BLUE:
                return (-mmt[cube.r1][cube.g1][cube.b0]
                    + mmt[cube.r1][cube.g0][cube.b0]
                    + mmt[cube.r0][cube.g1][cube.b0]
                    - mmt[cube.r0][cube.g0][cube.b0]);
            default:
                return 0;
        }
    }

    private static long top(final Box cube, final int color, final int pos, final long[][][] mmt) {
        /* Compute remainder of Vol(cube, mmt), substituting pos for */
        /* r1, g1, or b1 (depending on color) */
        switch (color) {
            case RED:
                return (mmt[pos][cube.g1][cube.b1]
                    - mmt[pos][cube.g1][cube.b0]
                    - mmt[pos][cube.g0][cube.b1]
                    + mmt[pos][cube.g0][cube.b0]);
            case GREEN:
                return (mmt[cube.r1][pos][cube.b1]
                    - mmt[cube.r1][pos][cube.b0]
                    - mmt[cube.r0][pos][cube.b1]
                    + mmt[cube.r0][pos][cube.b0]);
            case BLUE:
                return (mmt[cube.r1][cube.g1][pos]
                    - mmt[cube.r1][cube.g0][pos]
                    - mmt[cube.r0][cube.g1][pos]
                    + mmt[cube.r0][cube.g0][pos]);
            default:
                return 0;
        }
    }

    private float variance(final Box cube) {
        /* Compute the weighted variance of a box */
        /* NB: as with the raw statistics, this is really the variance * size */
        final long dr = volume(cube, mr);
        final long dg = volume(cube, mg);
        final long db = volume(cube, mb);
        final float xx = m2[cube.r1][cube.g1][cube.b1]
            - m2[cube.r1][cube.g1][cube.b0]
            - m2[cube.r1][cube.g0][cube.b1]
            + m2[cube.r1][cube.g0][cube.b0]
            - m2[cube.r0][cube.g1][cube.b1]
            + m2[cube.r0][cube.g1][cube.b0]
            + m2[cube.r0][cube.g0][cube.b1]
            - m2[cube.r0][cube.g0][cube.b0];
        return xx - (dr * dr + dg * dg + db * db) / (float) volume(cube, wt);
    }

    /* We want to minimize the sum of the variances of two subboxes.
     * The sum(c^2) terms can be ignored since their sum over both subboxes
     * is the same (the sum for the whole box) no matter where we split.
     * The remaining terms have a minus sign in the variance formula,
     * so we drop the minus sign and MAXIMIZE the sum of the two terms.
     */
    private float maximize(final Box cube, final int color, final int first, final int last, final int[] cut,
                           final long wholeR, final long wholeG, final long wholeB, final long wholeW) {
        final long baseR = bottom(cube, color, mr);
        final long baseG = bottom(cube, color, mg);
        final long baseB = bottom(cube, color, mb);
        final long baseW = bottom(cube, color, wt);

        float max = 0.0f;
        cut[0] = -1;

        for (int i = first; i < last; ++i) {
            long halfR = baseR + top(cube, color, i, mr);
            long halfG = baseG + top(cube, color, i, mg);
            long halfB = baseB + top(cube, color, i, mb);
            long halfW = baseW + top(cube, color, i, wt);
            /* now half_x is sum over lower half of box, if split at i */
            if (halfW == 0) { /* subbox could be empty of pixels! */
                continue;    /* never split into an empty box */
            }
            float temp = (halfR * halfR + halfG * halfG + halfB * halfB) / (float) halfW;
            halfR = wholeR - halfR;
            halfG = wholeG - halfG;
            halfB = wholeB - halfB;
            halfW = wholeW - halfW;
            if (halfW == 0) { /* subbox could be empty of pixels! */
                continue; /* never split into an empty box */
            }
            temp += (halfR * halfR + halfG * halfG + halfB * halfB) / (float) halfW;

            if (temp > max) {
                max = temp;
                cut[0] = i;
            }
        }

        return max;
    }

    private boolean cut(final Box set1, final Box set2) {
        final int[] cutR = new int[1];
        final int[] cutG = new int[1];
        final int[] cutB = new int[1];

        final long wholeR = volume(set1, mr);
        final long wholeG = volume(set1, mg);
        final long wholeB = volume(set1, mb);
        final long wholeW = volume(set1, wt);

        final float maxr = maximize(set1, RED, set1.r0 + 1, set1.r1, cutR,
            wholeR, wholeG, wholeB, wholeW);
        final float maxg = maximize(set1, GREEN, set1.g0 + 1, set1.g1, cutG,
            wholeR, wholeG, wholeB, wholeW);
        final float maxb = maximize(set1, BLUE, set1.b0 + 1, set1.b1, cutB,
            wholeR, wholeG, wholeB, wholeW);

        final int color;
        if (maxr >= maxg && maxr >= maxb) {
            color = RED;
            if (cutR[0] < 0) {
                return false; /* can't split the box */
            }
        }
        else if (maxg >= maxr && maxg >= maxb) {
            color = GREEN;
        }
        else {
            color = BLUE;
        }
        set2.r1 = set1.r1;
        set2.g1 = set1.g1;
        set2.b1 = set1.b1;

        switch (color) {
            case RED:
                set2.r0 = set1.r1 = cutR[0];
                set2.g0 = set1.g0;
                set2.b0 = set1.b0;
                break;
            case GREEN:
                set2.g0 = set1.g1 = cutG[0];
                set2.r0 = set1.r0;
                set2.b0 = set1.b0;
                break;
            case BLUE:
                set2.b0 = set1.b1 = cutB[0];
                set2.r0 = set1.r0;
                set2.g0 = set1.g0;
                break;
        }
        set1.vol = (set1.r1 - set1.r0) * (set1.g1 - set1.g0) * (set1.b1 - set1.b0);
        set2.vol = (set2.r1 - set2.r0) * (set2.g1 - set2.g0) * (set2.b1 - set2.b0);

        return true;
    }

    private void mark(final Box cube, final int label, final int[] tag) {
        for (int r = cube.r0 + 1; r <= cube.r1; ++r) {
            for (int g = cube.g0 + 1; g <= cube.g1; ++g) {
                for (int b = cube.b0 + 1; b <= cube.b1; ++b) {
                    tag[(r << 10) + (r << 6) + r + (g << 5) + g + b] = label;
                }
            }
        }
    }

    private static final class Box {
        int r0;     /* min value, exclusive */
        int r1;     /* max value, inclusive */
        int g0;
        int g1;
        int b0;
        int b1;
        int vol;
    }
}
