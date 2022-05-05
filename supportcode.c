ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];

ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];