ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);

ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, 0, k, a)] + 9.0 * ff[fineindex(i, 2, k, a)] - ff[fineindex(i, 4, k, a)]) / 16.0;

ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(nx, j, k, a)] + 9.0 * ff[fineindex(nx - 2, j, k, a)] - ff[fineindex(nx - 4, j, k, a)]) / 16.0;

ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, ny, k, a)] + 9.0 * ff[fineindex(i, ny - 2, k, a)] - ff[fineindex(i, ny - 4, k, a)]) / 16.0;

ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(0, j, k, a)] + 9.0 * ff[fineindex(2, j, k, a)] - ff[fineindex(4, j, k, a)]) / 16.0;

ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, j, 0, a)] + 9.0 * ff[fineindex(i, j, 2, a)] - ff[fineindex(i, j, 4, a)]) / 16.0;

ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, j, nz, a)] + 9.0 * ff[fineindex(i, j, nz - 2, a)] - ff[fineindex(i, j, nz - 4, a)]) / 16.0;

ff[fineindex(i, j, k, a)] = (ff[fineindex(i - 3, j, k + 3, a)] - 9.0 * ff[fineindex(i - 3, j, k + 1, a)] - 9.0 * ff[fineindex(i - 3, j, k - 1, a)] + ff[fineindex(i - 3, j, k - 3, a)] - 9.0 * ff[fineindex(i - 1, j, k + 3, a)] + 81.0 * ff[fineindex(i - 1, j, k + 1, a)] + 81.0 * ff[fineindex(i - 1, j, k - 1, a)] - 9.0 * ff[fineindex(i - 1, j, k - 3, a)] - 9.0 * ff[fineindex(i + 1, j, k + 3, a)] + 81.0 * ff[fineindex(i + 1, j, k + 1, a)] + 81.0 * ff[fineindex(i + 1, j, k - 1, a)] - 9.0 * ff[fineindex(i + 1, j, k - 3, a)] + ff[fineindex(i + 3, j, k + 3, a)] - 9.0 * ff[fineindex(i + 3, j, k + 1, a)] - 9.0 * ff[fineindex(i + 3, j, k - 1, a)] + ff[fineindex(i + 3, j, k - 3, a)]) / 256.0;

ff[fineindex(i, j, k, a)] = (ff[fineindex(i, j - 3, k + 3, a)] - 9.0 * ff[fineindex(i, j - 3, k + 1, a)] - 9.0 * ff[fineindex(i, j - 3, k - 1, a)] + ff[fineindex(i, j - 3, k - 3, a)] - 9.0 * ff[fineindex(i, j - 1, k + 3, a)] + 81.0 * ff[fineindex(i, j - 1, k + 1, a)] + 81.0 * ff[fineindex(i, j - 1, k - 1, a)] - 9.0 * ff[fineindex(i, j - 1, k - 3, a)] - 9.0 * ff[fineindex(i, j + 1, k + 3, a)] + 81.0 * ff[fineindex(i, j + 1, k + 1, a)] + 81.0 * ff[fineindex(i, j + 1, k - 1, a)] - 9.0 * ff[fineindex(i, j + 1, k - 3, a)] + ff[fineindex(i, j + 3, k + 3, a)] - 9.0 * ff[fineindex(i, j + 3, k + 1, a)] - 9.0 * ff[fineindex(i, j + 3, k - 1, a)] + ff[fineindex(i, j + 3, k - 3, a)]) / 256.0;

ff[fineindex(i, j, k, a)] = (ff[fineindex(i - 3, j, k + 3, a)] - 9.0 * ff[fineindex(i - 3, j, k + 1, a)] - 9.0 * ff[fineindex(i - 3, j, k - 1, a)] + ff[fineindex(i - 3, j, k - 3, a)] - 9.0 * ff[fineindex(i - 1, j, k + 3, a)] + 81.0 * ff[fineindex(i - 1, j, k + 1, a)] + 81.0 * ff[fineindex(i - 1, j, k - 1, a)] - 9.0 * ff[fineindex(i - 1, j, k - 3, a)] - 9.0 * ff[fineindex(i + 1, j, k + 3, a)] + 81.0 * ff[fineindex(i + 1, j, k + 1, a)] + 81.0 * ff[fineindex(i + 1, j, k - 1, a)] - 9.0 * ff[fineindex(i + 1, j, k - 3, a)] + ff[fineindex(i + 3, j, k + 3, a)] - 9.0 * ff[fineindex(i + 3, j, k + 1, a)] - 9.0 * ff[fineindex(i + 3, j, k - 1, a)] + ff[fineindex(i + 3, j, k - 3, a)]) / 256.0;

ff[fineindex(i, j, k, a)] = (ff[fineindex(i, j - 3, k + 3, a)] - 9.0 * ff[fineindex(i, j - 3, k + 1, a)] - 9.0 * ff[fineindex(i, j - 3, k - 1, a)] + ff[fineindex(i, j - 3, k - 3, a)] - 9.0 * ff[fineindex(i, j - 1, k + 3, a)] + 81.0 * ff[fineindex(i, j - 1, k + 1, a)] + 81.0 * ff[fineindex(i, j - 1, k - 1, a)] - 9.0 * ff[fineindex(i, j - 1, k - 3, a)] - 9.0 * ff[fineindex(i, j + 1, k + 3, a)] + 81.0 * ff[fineindex(i, j + 1, k + 1, a)] + 81.0 * ff[fineindex(i, j + 1, k - 1, a)] - 9.0 * ff[fineindex(i, j + 1, k - 3, a)] + ff[fineindex(i, j + 3, k + 3, a)] - 9.0 * ff[fineindex(i, j + 3, k + 1, a)] - 9.0 * ff[fineindex(i, j + 3, k - 1, a)] + ff[fineindex(i, j + 3, k - 3, a)]) / 256.0;

ff[fineindex(i, j, k, a)] = (ff[fineindex(i - 3, j + 3, k, a)] - 9.0 * ff[fineindex(i - 3, j + 1, k, a)] - 9.0 * ff[fineindex(i - 3, j - 1, k, a)] + ff[fineindex(i - 3, j - 3, k, a)] - 9.0 * ff[fineindex(i - 1, j + 3, k, a)] + 81.0 * ff[fineindex(i - 1, j + 1, k, a)] + 81.0 * ff[fineindex(i - 1, j - 1, k, a)] - 9.0 * ff[fineindex(i - 1, j - 3, k, a)] - 9.0 * ff[fineindex(i + 1, j + 3, k, a)] + 81.0 * ff[fineindex(i + 1, j + 1, k, a)] + 81.0 * ff[fineindex(i + 1, j - 1, k, a)] - 9.0 * ff[fineindex(i + 1, j - 3, k, a)] + ff[fineindex(i + 3, j + 3, k, a)] - 9.0 * ff[fineindex(i + 3, j + 1, k, a)] - 9.0 * ff[fineindex(i + 3, j - 1, k, a)] + ff[fineindex(i + 3, j - 3, k, a)]) / 256.0;

ff[fineindex(i, j, k, a)] = (ff[fineindex(i - 3, j + 3, k, a)] - 9.0 * ff[fineindex(i - 3, j + 1, k, a)] - 9.0 * ff[fineindex(i - 3, j - 1, k, a)] + ff[fineindex(i - 3, j - 3, k, a)] - 9.0 * ff[fineindex(i - 1, j + 3, k, a)] + 81.0 * ff[fineindex(i - 1, j + 1, k, a)] + 81.0 * ff[fineindex(i - 1, j - 1, k, a)] - 9.0 * ff[fineindex(i - 1, j - 3, k, a)] - 9.0 * ff[fineindex(i + 1, j + 3, k, a)] + 81.0 * ff[fineindex(i + 1, j + 1, k, a)] + 81.0 * ff[fineindex(i + 1, j - 1, k, a)] - 9.0 * ff[fineindex(i + 1, j - 3, k, a)] + ff[fineindex(i + 3, j + 3, k, a)] - 9.0 * ff[fineindex(i + 3, j + 1, k, a)] - 9.0 * ff[fineindex(i + 3, j - 1, k, a)] + ff[fineindex(i + 3, j - 3, k, a)]) / 256.0;