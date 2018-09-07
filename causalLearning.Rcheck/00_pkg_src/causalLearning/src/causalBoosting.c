#include <math.h>       // for pow()
#include <Rinternals.h> // for NA_REAL
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  // for timing

const double eps = 1e-8;

double mean(double *x, int length) {
  double sum;
  int i;
  sum = 0;
  for (i = 0; i < length; i++) {
    sum += x[i];
  }
  return(sum/length);
}

double meanWtd(double *x, int length, double *w) {
  double den, sum;
  int i;
  den = 0;
  sum = 0;
  for (i = 0; i < length; i++) {
    den += w[i];
    sum += w[i]*x[i];
  }
  return(sum/den);
}

double sampVar(double *x, int length, double mean) {
  double sumsq;
  int i;
  sumsq = 0;
  for (i = 0; i < length; i++) {
    sumsq += pow((x[i] - mean), 2);
  }
  return(sumsq/(length-1));
}

double sampVarWtd(double *x, int length, double mean, double *w) {
  double den, sumsq;
  int i;
  den = 0;
  sumsq = 0;
  for (i = 0; i < length; i++) {
    den += w[i];
    sumsq += w[i]*pow((x[i] - mean), 2);
  }
  return(sumsq/(den-1));
}

int meanWtdStrat(double *x, int length, double *w, int *s, int ns, double *totWeight, double *mean) {
  // s stores the index of the bin that each x belongs to
  for (int k = 0; k < ns; ++k) {
    mean[k] = 0;
    totWeight[k] = 0;
  }
  for (int i = 0; i < length; ++i) {
    mean[s[i]] += w[i] * x[i];
    totWeight[s[i]] += w[i];
  }
  for (int k = 0; k < ns; ++k) {
    if (totWeight[k] > 0)
      mean[k] /= totWeight[k];  
    else
      mean[k] = 0;
  }
  return 1; 
}

int sampVarWtdStrat(double *x, int length, double *mean, double *w, double *totWeight, int *s, int ns, double *var) {
  // s stores the index of the bin that each x belongs to
  for (int k = 0; k < ns; ++k) {
    var[k] = 0;
  }
  for (int i = 0; i < length; ++i) {
    var[s[i]] += w[i] * (x[i] - mean[s[i]]) * (x[i] - mean[s[i]]);
  }
  for (int k = 0; k < ns; ++k) {
    if (totWeight[k] > 1)
      var[k] /= (totWeight[k] - 1);
    else 
      var[k] = 0;
  }
  return 1; 
}

int constVarStrat(double *x, int length, double *mean, double *w, double *totWeight, int *s, int ns, double *var) {
  // used when we don't want to estimate sample variance of the treated and untreated
  // since we are computing the relative splitting criterion, actual number set doesn't matter
  // and we juse set it to 1.
  for (int k = 0; k < ns; ++k) {
    var[k] = 1;
  }
  return 1;
}

double splitCriterion(double *y, int *tx, double *w, int *a, int *b, int *n)
{
  double txAmean, cxAmean, txBmean, cxBmean;
  double txAvar, cxAvar, txBvar, cxBvar;
  double tauA, tauB, varA, varB;
  int nn = *n;
  double *txA = (double *)malloc(nn * sizeof(double));
  double *cxA = (double *)malloc(nn * sizeof(double));
  double *txB = (double *)malloc(nn * sizeof(double));
  double *cxB = (double *)malloc(nn * sizeof(double));
  double *txAw = (double *)malloc(nn * sizeof(double));
  double *cxAw = (double *)malloc(nn * sizeof(double));
  double *txBw = (double *)malloc(nn * sizeof(double));
  double *cxBw = (double *)malloc(nn * sizeof(double));
  int txAlen, cxAlen, txBlen, cxBlen, i;
  txAlen = 0;
  cxAlen = 0;
  txBlen = 0;
  cxBlen = 0;
  for (i = 0; i < *n; i++) {
    if (tx[i] && a[i]) {
      txA[txAlen] = y[i];
      txAw[txAlen] = w[i];
      txAlen += 1;
    }
    else if (tx[i] == 0 && a[i]) {
      cxA[cxAlen] = y[i];
      cxAw[cxAlen] = w[i];
      cxAlen += 1;
    }
    else if (tx[i] && b[i]) {
      txB[txBlen] = y[i];
      txBw[txBlen] = w[i];
      txBlen += 1;
    }
    else if (tx[i] == 0 && b[i]) {
      cxB[cxBlen] = y[i];
      cxBw[cxBlen] = w[i];
      cxBlen += 1;
    }
  }
//  if ((txAlen + cxAlen < 5) || (txBlen + cxBlen < 5)) {
//    return(0);
  //  }

  double out;

  if (txAlen < 2 || cxAlen < 2 || txBlen < 2 || cxBlen < 2) {
    out = 0;
  } else {
    txAmean = meanWtd(txA, txAlen, txAw);
    cxAmean = meanWtd(cxA, cxAlen, cxAw);
    txBmean = meanWtd(txB, txBlen, txBw);
    cxBmean = meanWtd(cxB, cxBlen, cxBw);
    txAvar = sampVarWtd(txA, txAlen, txAmean, txAw);
    cxAvar = sampVarWtd(cxA, cxAlen, cxAmean, cxAw);
    txBvar = sampVarWtd(txB, txBlen, txBmean, txBw);
    cxBvar = sampVarWtd(cxB, cxBlen, cxBmean, cxBw);
    tauA = txAmean - cxAmean;
    tauB = txBmean - cxBmean;
    varA = (txAvar / txAlen) + (cxAvar/ cxAlen);
    varB = (txBvar / txBlen) + (cxBvar/ cxBlen);
    out = fabs(tauA - tauB) / sqrt(varA + varB);
  }

  free(txA);
  free(cxA);
  free(txB);
  free(cxB);
  free(txAw);
  free(cxAw);
  free(txBw);
  free(cxBw);

  return out;
}


double splitCriterion_propensity(double *y, int *tx, double *w, int *a, int *b, int *n,
                                 int *s, int *ns, int *isConstVar)
// s should be 0-based
{
  double *txAmean = (double *)malloc((*ns) * sizeof(double));
  double *cxAmean = (double *)malloc((*ns) * sizeof(double));
  double *txBmean = (double *)malloc((*ns) * sizeof(double));
  double *cxBmean = (double *)malloc((*ns) * sizeof(double));
  double *txAvar = (double *)malloc((*ns) * sizeof(double));
  double *cxAvar = (double *)malloc((*ns) * sizeof(double));
  double *txBvar = (double *)malloc((*ns) * sizeof(double));
  double *cxBvar = (double *)malloc((*ns) * sizeof(double));
  double *txAtotw = (double *)malloc((*ns) * sizeof(double));
  double *cxAtotw = (double *)malloc((*ns) * sizeof(double));
  double *txBtotw = (double *)malloc((*ns) * sizeof(double));
  double *cxBtotw = (double *)malloc((*ns) * sizeof(double));
  double tauA = 0, tauB = 0, varA = 0, varB = 0;
  double twA = 0, twB = 0;

  int nn = *n;
  double *txA = (double *)malloc(nn * sizeof(double));
  double *cxA = (double *)malloc(nn * sizeof(double));
  double *txB = (double *)malloc(nn * sizeof(double));
  double *cxB = (double *)malloc(nn * sizeof(double));
  double *txAw = (double *)malloc(nn * sizeof(double));
  double *cxAw = (double *)malloc(nn * sizeof(double));
  double *txBw = (double *)malloc(nn * sizeof(double));
  double *cxBw = (double *)malloc(nn * sizeof(double));
  int *txAs = (int *)malloc(nn * sizeof(int));
  int *cxAs = (int *)malloc(nn * sizeof(int));
  int *txBs = (int *)malloc(nn * sizeof(int));
  int *cxBs = (int *)malloc(nn * sizeof(int));

  int txAlen = 0, cxAlen = 0, txBlen = 0, cxBlen = 0;
  
  for (int i = 0; i < *n; i++) {
    if (tx[i] && a[i]) {
      txA[txAlen] = y[i];
      txAw[txAlen] = w[i];
      txAs[txAlen] = s[i];
      txAlen += 1;
    }
    else if (tx[i] == 0 && a[i]) {
      cxA[cxAlen] = y[i];
      cxAw[cxAlen] = w[i];
      cxAs[cxAlen] = s[i];
      cxAlen += 1;
    }
    else if (tx[i] && b[i]) {
      txB[txBlen] = y[i];
      txBw[txBlen] = w[i];
      txBs[txBlen] = s[i];
      txBlen += 1;
    }
    else if (tx[i] == 0 && b[i]) {
      cxB[cxBlen] = y[i];
      cxBw[cxBlen] = w[i];
      cxBs[cxBlen] = s[i];
      cxBlen += 1;
    }
  }
  
  meanWtdStrat(txA, txAlen, txAw, txAs, *ns, txAtotw, txAmean);
  meanWtdStrat(cxA, cxAlen, cxAw, cxAs, *ns, cxAtotw, cxAmean);
  meanWtdStrat(txB, txBlen, txBw, txBs, *ns, txBtotw, txBmean);
  meanWtdStrat(cxB, cxBlen, cxBw, cxBs, *ns, cxBtotw, cxBmean);
  
  if (!(*isConstVar)) {
    sampVarWtdStrat(txA, txAlen, txAmean, txAw, txAtotw, txAs, *ns, txAvar);
    sampVarWtdStrat(cxA, cxAlen, cxAmean, cxAw, cxAtotw, cxAs, *ns, cxAvar);
    sampVarWtdStrat(txB, txBlen, txBmean, txBw, txBtotw, txBs, *ns, txBvar);
    sampVarWtdStrat(cxB, cxBlen, cxBmean, cxBw, cxBtotw, cxBs, *ns, cxBvar);
  }
  else {
    constVarStrat(txA, txAlen, txAmean, txAw, txAtotw, txAs, *ns, txAvar);
    constVarStrat(cxA, cxAlen, cxAmean, cxAw, cxAtotw, cxAs, *ns, cxAvar);
    constVarStrat(txB, txBlen, txBmean, txBw, txBtotw, txBs, *ns, txBvar);
    constVarStrat(cxB, cxBlen, cxBmean, cxBw, cxBtotw, cxBs, *ns, cxBvar);
  }
  
  for (int k = 0; k < *ns; ++k) {
    int wsA = txAtotw[k] + cxAtotw[k], wsB = txBtotw[k] + cxBtotw[k];
    if (txAtotw[k] > eps && cxAtotw[k] > eps) {
      if ((*isConstVar) || (txAtotw[k] > 1 && cxAtotw[k] > 1)) {
        tauA += wsA * (txAmean[k] - cxAmean[k]);
        twA += wsA;
        varA += wsA * wsA * (txAvar[k]/txAtotw[k] + cxAvar[k]/cxAtotw[k]);
      }
    }
    if (txBtotw[k] > eps && cxBtotw[k] > eps) {
      if ((*isConstVar) || (txBtotw[k] > 1 && cxBtotw[k] > 1)) {
        tauB += wsB * (txBmean[k] - cxBmean[k]);
        twB += wsB;
        varB += wsB * wsB * (txBvar[k]/txBtotw[k] + cxBvar[k]/cxBtotw[k]);
      }
    }
  }

  free(txAmean);
  free(cxAmean);
  free(txBmean);
  free(cxBmean);
  free(txAvar);
  free(cxAvar);
  free(txBvar);
  free(cxBvar);
  free(txAtotw);
  free(cxAtotw);
  free(txBtotw);
  free(cxBtotw);

  free(txA);
  free(cxA);
  free(txB);
  free(txB);
  free(txAw);
  free(cxAw);
  free(txBw);
  free(cxBw);
  free(txAs);
  free(cxAs);
  free(txBs);
  free(cxBs);
  
  if (fabs(varA) < eps || fabs(varB) < eps) return 0;

  if (twA != 0) {
    tauA /= twA;
    varA /= twA * twA;
  }
  else {
    return 0;
  }
  
  if (twB != 0) {
    tauB /= twB;
    varB /= twB * twB;
  }
  else {
    return 0;
  }

  return (fabs(tauA - tauB) / sqrt(varA + varB));


  // txAmean = meanWtd(txA, txAlen, txAw);
  // cxAmean = meanWtd(cxA, cxAlen, cxAw);
  // txBmean = meanWtd(txB, txBlen, txBw);
  // cxBmean = meanWtd(cxB, cxBlen, cxBw);
  // txAvar = sampVarWtd(txA, txAlen, txAmean, txAw);
  // cxAvar = sampVarWtd(cxA, cxAlen, cxAmean, cxAw);
  // txBvar = sampVarWtd(txB, txBlen, txBmean, txBw);
  // cxBvar = sampVarWtd(cxB, cxBlen, cxBmean, cxBw);
  // tauA = txAmean - cxAmean;
  // tauB = txBmean - cxBmean;
  // varA = (txAvar / txAlen) + (cxAvar/ cxAlen);
  // varB = (txBvar / txBlen) + (cxBvar/ cxBlen);
  // return(fabs(tauA - tauB) / sqrt(varA + varB));
}


void findSplit(double *X, double *y, int *tx, double *w, int *order,
  int *n, int *q, int *variable, double *value, double *objective)
{
  int i, j, above[*n], below[*n], variableMax;
  double valueMax, criterionMax, criterion, previous;
  criterionMax = 0;
  for (j = 0; j < *q; j++) {
    for (i = 0; i < *n; i++) {
      above[i] = 0;
      if (w[i] > 0) {
        above[i] = 1;
      }
      below[i] = 0;
    }
    previous = X[order[j * *n] + j * *n];
    for (i = 0; i < *n-1; i++) {
      if (w[order[i + j * *n]] > 0) {
        if (X[order[i + j * *n] + j * *n] != previous) {
          criterion = splitCriterion(y, tx, w, above, below, n);
          if (criterion > criterionMax) {
            criterionMax = criterion;
            variableMax = j;
            valueMax = X[order[i + j * *n] + j * *n];
          }
        }
        above[order[i + j * *n]] = 0;
        below[order[i + j * *n]] = 1;
        previous = X[order[i + j * *n] + j * *n];
      }
    }
  }
  *variable = variableMax;
  *value = valueMax;
  *objective = criterionMax;
  if (criterionMax == 0) {
    *variable = NA_REAL;
    *value = NA_REAL;
    *objective = NA_REAL;
  }
}

int comp_int(const void *a,const void *b) {
  int *x = (int *) a;
  int *y = (int *) b;
  return *x - *y;
}

// find unique numbers in first n elements of u -> saved in the first k slots, return k
int findUnique(int *u, int n) {
  if (n == 0) return 0;
  qsort(u, n, sizeof(int), comp_int);
  int k = 1, prev = u[0];
  for (int i = 1; i < n; ++i) {
    if (u[i] != prev) {
      u[k] = u[i];
      prev = u[i];
      k++;
    }
  }
  return k;
}


// consider every splitSpread * nn points to the left node
void findBestSplit(double *X, double *y, int *tx, double *w, int *order,
 int *n, int *q, int *leaf, int *leaf_var, double *leaf_val, double *leaf_criterion,
 int *node, int *variable, double *value, double *objective, double *splitSpread) {

  int nn = *n;

  int *uniqueLeaf = (int *)malloc(nn * sizeof(int));
  // int uniqueLeaf[nn];
  memcpy(uniqueLeaf, leaf, nn * sizeof(int));
  int numLeaves = findUnique(uniqueLeaf, nn);
  
  int variableMax[numLeaves];
  double valueMax[numLeaves];
  double criterionMax[numLeaves];
  
  for (int k = 0; k < numLeaves; ++k) {
    int currLeaf = uniqueLeaf[k];
    if (leaf_criterion[currLeaf] > -eps) {  // only applicable when splitCriterion >= 0
      criterionMax[k] = leaf_criterion[currLeaf];
      variableMax[k] = leaf_var[currLeaf];
      valueMax[k] = leaf_val[currLeaf];
    }
    else {
      double *wk = (double *)malloc(nn * sizeof(double));
      // double wk[nn];
      int currSize = 0;
      memcpy(wk, w, nn * sizeof(double));
      for (int i = 0; i < nn; ++i) {
        if (leaf[i] != currLeaf) {
          wk[i] = 0;
        }
        else {
          currSize++;
        }
      }

      int splitGap = (int)round(currSize * (*splitSpread));
      if (splitGap < 1) {
        splitGap = 1;
      }

      criterionMax[k] = 0;
      variableMax[k] = 0;
      for (int j = 0; j < *q; j++) {
        int above[nn], below[nn];
        for (int i = 0; i < nn; i++) {
          above[i] = 0;
          if (wk[i] > eps) {
            above[i] = 1;
          }
          below[i] = 0;
        }

        // double previous = X[order[j * *n] + j * *n] - 1;
        // for (int i = 0; i < *n-1; i++) {
        //   if (wk[order[i + j * *n]] > eps) {
        //     if (fabs(X[order[i + j * *n] + j * *n] - previous) > eps) {
        //       double criterion = splitCriterion(y, tx, wk, above, below, n);
        //       if (criterion > criterionMax[k]) {
        //         criterionMax[k] = criterion;
        //         variableMax[k] = j;
        //         valueMax[k] = X[order[i + j * *n] + j * *n];
        //       }
        //     }
        //     above[order[i + j * *n]] = 0;
        //     below[order[i + j * *n]] = 1;
        //     previous = X[order[i + j * *n] + j * *n];
        //   }
        // }
        double current;
        int splitCount = 0, kks = -1;
        for (int i = 0; i < nn; i++) {
          if (wk[order[i + j * nn]] > eps) {
            // if (fabs(X[order[i + j * nn] + j * nn] - previous) > eps) {
            if (splitCount / splitGap > kks) {
              double criterion = splitCriterion(y, tx, wk, above, below, n);
              kks = splitCount / splitGap;
              if (criterion > criterionMax[k]) {
                criterionMax[k] = criterion;
                variableMax[k] = j;
                valueMax[k] = X[order[i + j * nn] + j * nn];
              }
            }
            // }
            current = X[order[i + j * nn] + j * nn];
            for (int u = i; u < nn; ++u) {
              if (wk[order[u + j * nn]] > eps) {
                if (fabs(X[order[u + j * nn] + j * nn] - current) <= eps) {
                  splitCount++;
                  above[order[u + j * nn]] = 0;
                  below[order[u + j * nn]] = 1;
                }
                else {
                  i = u - 1;
                  break;
                }
              }
            }
            // above[order[i + j * nn]] = 0;
            // below[order[i + j * nn]] = 1;
          }
        }
      }
      leaf_criterion[currLeaf] = criterionMax[k];
      leaf_var[currLeaf] = variableMax[k];
      leaf_val[currLeaf] = valueMax[k];
      free(wk);
    }
  }

  double optCriterion = 0;  // only when splitCriterion >= 0
  int optIdx = 0;
  for (int k = 0; k < numLeaves; ++k) {
    if (criterionMax[k] > optCriterion) {
      optIdx = k;
      optCriterion = criterionMax[k];
    }
  }

  *variable = variableMax[optIdx];
  *value = valueMax[optIdx];
  *objective = criterionMax[optIdx];
  *node = uniqueLeaf[optIdx];

  if (optCriterion == 0) {
   *variable = NA_INTEGER;
   *value = NA_REAL;
   *objective = NA_REAL;
   *node = NA_INTEGER;
  }

  free(uniqueLeaf);
  // end = clock();
  // *cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}


// void findBestSplit_debug(double *X, double *y, int *tx, double *w, int *order,
//  int *n, int *q, int *leaf, int *leaf_var, double *leaf_val, double *leaf_criterion,
//  int *node, int *variable, double *value, double *objective, double *cpu_time_used, int *count)
// {
//   clock_t start, end;

//   // start = clock();

//   int uniqueLeaf[*n];
//   memcpy(uniqueLeaf, leaf, (*n) * sizeof(int));
//   int numLeaves = findUnique(uniqueLeaf, *n);
  
//   int variableMax[numLeaves];
//   double valueMax[numLeaves];
//   double criterionMax[numLeaves];
  
//   for (int k = 0; k < numLeaves; ++k) {
//     int currLeaf = uniqueLeaf[k];
//     if (leaf_criterion[currLeaf] > -eps) {  // only applicable when splitCriterion >= 0
//       criterionMax[k] = leaf_criterion[currLeaf];
//       variableMax[k] = leaf_var[currLeaf];
//       valueMax[k] = leaf_val[currLeaf];
//     }
//     else {
//       double wk[*n];
//       memcpy(wk, w, (*n) * sizeof(double));
//       for (int i = 0; i < *n; ++i) {
//         if (leaf[i] != currLeaf) {
//           wk[i] = 0;
//         }
//       }

//       criterionMax[k] = 0;
//       variableMax[k] = 0;
//       for (int j = 0; j < *q; j++) {
//         int above[*n], below[*n];
//         double criterion, previous;
//         for (int i = 0; i < *n; i++) {
//           above[i] = 0;
//           if (wk[i] > eps) {
//             above[i] = 1;
//           }
//           below[i] = 0;
//         }
//         previous = X[order[j * *n] + j * *n] - 1;
//         for (int i = 0; i < *n-1; i++) {
//           if (wk[order[i + j * *n]] > eps) {
//             if (fabs(X[order[i + j * *n] + j * *n] - previous) > eps) {

//       *count = (*count) + 1;
//       start = clock();
//               criterion = splitCriterion(y, tx, wk, above, below, n);
//     end = clock();
//     *cpu_time_used += ((double) (end - start)) / CLOCKS_PER_SEC;
//               if (criterion > criterionMax[k]) {
//                 criterionMax[k] = criterion;
//                 variableMax[k] = j;
//                 valueMax[k] = X[order[i + j * *n] + j * *n];
//               }
//             }
//             above[order[i + j * *n]] = 0;
//             below[order[i + j * *n]] = 1;
//             previous = X[order[i + j * *n] + j * *n];
//           }
//         }
//       }
//       leaf_criterion[currLeaf] = criterionMax[k];
//       leaf_var[currLeaf] = variableMax[k];
//       leaf_val[currLeaf] = valueMax[k];
//     }
//   }

//   double optCriterion = 0;  // only when splitCriterion >= 0
//   int optIdx = 0;
//   for (int k = 0; k < numLeaves; ++k) {
//     if (criterionMax[k] > optCriterion) {
//       optIdx = k;
//       optCriterion = criterionMax[k];
//     }
//   }

//   *variable = variableMax[optIdx];
//   *value = valueMax[optIdx];
//   *objective = criterionMax[optIdx];
//   *node = uniqueLeaf[optIdx];

//   if (optCriterion == 0) {
//    *variable = NA_INTEGER;
//    *value = NA_REAL;
//    *objective = NA_REAL;
//    *node = NA_INTEGER;
//   }

//   // end = clock();
//   // *cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
// }


void findBestSplit_propensity(double *X, double *y, int *tx, double *w,
  int *order,
 int *n, int *q, int *leaf, int *leaf_var, double *leaf_val, double *leaf_criterion, int *s, int *ns, int *isConstVar,
 int *node, int *variable, double *value, double *objective, double *splitSpread)
{
  int nn = *n;
  int *uniqueLeaf = (int *)malloc(nn * sizeof(int));
  memcpy(uniqueLeaf, leaf, nn * sizeof(int));
  int numLeaves = findUnique(uniqueLeaf, nn);
  
  int variableMax[numLeaves];
  double valueMax[numLeaves];
  double criterionMax[numLeaves];
  
  for (int k = 0; k < numLeaves; ++k) {
    int currLeaf = uniqueLeaf[k];
    if (leaf_criterion[currLeaf] > -eps) {  // only applicable when splitCriterion >= 0
      criterionMax[k] = leaf_criterion[currLeaf];
      variableMax[k] = leaf_var[currLeaf];
      valueMax[k] = leaf_val[currLeaf];
    }
    else {
      double *wk = (double *)malloc(nn * sizeof(double));
      memcpy(wk, w, nn * sizeof(double));
      int currSize = 0;
      for (int i = 0; i < nn; ++i) {
        if (leaf[i] != currLeaf) {
          wk[i] = 0;
        }
        else {
          currSize++;
        }
      }

      int splitGap = (int)round(currSize * (*splitSpread));
      if (splitGap < 1) {
        splitGap = 1;
      }

      criterionMax[k] = 0;
      variableMax[k] = 0;
      for (int j = 0; j < *q; j++) {
        int above[nn], below[nn];
        for (int i = 0; i < nn; i++) {
          above[i] = 0;
          if (wk[i] > eps) {
            above[i] = 1;
          }
          below[i] = 0;
        }

        // double previous = X[order[j * *n] + j * *n] - 1;
        // for (int i = 0; i < *n-1; i++) {
        //   if (wk[order[i + j * *n]] > eps) {
        //     if (fabs(X[order[i + j * *n] + j * *n] - previous) > eps) {
        //       double criterion = splitCriterion_propensity(y, tx, wk, above, below, n, s, ns, isConstVar);
        //       if (criterion > criterionMax[k]) {
        //         criterionMax[k] = criterion;
        //         variableMax[k] = j;
        //         valueMax[k] = X[order[i + j * *n] + j * *n];
        //       }
        //     }
        //     above[order[i + j * *n]] = 0;
        //     below[order[i + j * *n]] = 1;
        //     previous = X[order[i + j * *n] + j * *n];
        //   }
        // }

        
        double current;
        int splitCount = 0, kks = -1;
        for (int i = 0; i < nn; i++) {
          if (wk[order[i + j * nn]] > eps) {
            // if (fabs(X[order[i + j * nn] + j * nn] - previous) > eps) {
            if (splitCount / splitGap > kks) {
              double criterion = splitCriterion_propensity(y, tx, wk, above, below, n, s, ns, isConstVar);
              kks = splitCount / splitGap;
              if (criterion > criterionMax[k]) {
                criterionMax[k] = criterion;
                variableMax[k] = j;
                valueMax[k] = X[order[i + j * nn] + j * nn];
              }
            }
            // }
            current = X[order[i + j * nn] + j * nn];
            for (int u = i; u < nn; ++u) {
              if (wk[order[u + j * nn]] > eps) {
                if (fabs(X[order[u + j * nn] + j * nn] - current) <= eps) {
                  splitCount++;
                  above[order[u + j * nn]] = 0;
                  below[order[u + j * nn]] = 1;
                }
                else {
                  i = u - 1;
                  break;
                }
              }
            }
            // above[order[i + j * nn]] = 0;
            // below[order[i + j * nn]] = 1;
          }
        }
      }
      leaf_criterion[currLeaf] = criterionMax[k];
      leaf_var[currLeaf] = variableMax[k];
      leaf_val[currLeaf] = valueMax[k];
      free(wk);
    }
  }

  double optCriterion = 0;  // only when splitCriterion >= 0
  int optIdx = -1;
  for (int k = 0; k < numLeaves; ++k) {
    if (criterionMax[k] > optCriterion) {
      optIdx = k;
      optCriterion = criterionMax[k];
    }
  }

  if (optIdx == -1) {
    *variable = NA_INTEGER;
    *value = NA_REAL;
    *objective = NA_REAL;
    *node = NA_INTEGER;
  } else {
    *variable = variableMax[optIdx];
    *value = valueMax[optIdx];
    *objective = criterionMax[optIdx];
    *node = uniqueLeaf[optIdx];
  }
  free(uniqueLeaf);
}

/* **********************************************************/
/* Prepare for causal boosting 4, main function under test. */
/* **********************************************************/




double meanByFactor(double *y, int *factor, int level, int size) {
  double sum = 0;
  int m = 0;
  for (int i = 0; i < size; ++i) {
    if (factor[i] == level) {
      sum += y[i];
      m++;
    }
  }
  if (m == 0) return 0;
  return sum / m;
}




int pmax_int(int *y, int size) {
  // size > 0
  int ymax = y[0];
  for (int i = 0; i < size; ++i) {
    if (y[i] > ymax) {
      ymax = y[i];
    }
  }
  return ymax;
}

struct num_ind {
  double val;
  double ind;
};

int cmp_num_ind(const void *a, const void *b) {
  struct num_ind *x = (struct num_ind *)a;
  struct num_ind *y = (struct num_ind *)b;
  return (x->val < y->val ? -1 : 1);
}

void order_index_double(double *y, int *order, int size) {
  // size > 0
  struct num_ind expanded[size];
  for (int i = 0; i < size; ++i) {
    expanded[i].val = y[i];
    expanded[i].ind = i;
  }
  qsort(expanded, size, sizeof(struct num_ind), cmp_num_ind);
  for (int i = 0; i < size; ++i) {
    order[i] = expanded[i].ind;
  }
  return;
}

void mean_propensity_tx(double *y, int *tx, double *w, int *leaf, int ell, int *n, int *s, int *ns, double *pred0, double *pred1) {
  double s_count[*ns][2];
  double ys_mean[*ns][2];

  for (int i = 0; i < *ns; ++i) {
    for (int j = 0; j < 2; ++j) {
      s_count[i][j] = 0;
      ys_mean[i][j] = 0;
    }
  }

  double count = 0;

  *pred0 = 0;
  *pred1 = 0;

  for (int j = 0; j < *ns; ++j) 
    for (int k = 0; k < 2; ++k)
      s_count[j][k] = 0;

  for (int i = 0; i < *n; ++i) {
    if (leaf[i] == ell) {
      s_count[s[i]][tx[i]] += w[i];
      ys_mean[s[i]][tx[i]] += w[i] * y[i];
    }
  }

  for (int j = 0; j < *ns; ++j)
    for (int k = 0; k < 2; ++k)
      if (s_count[j][k] > 0) {
        count += s_count[j][k];
        ys_mean[j][k] /= s_count[j][k];  // ns * 2 table entry-wise mean
      }

  for (int j = 0; j < *ns; ++j) {
    int nn = 0;
    for (int k = 0; k < 2; ++k) {
      nn += s_count[j][k];
    }
    *pred0 += nn * ys_mean[j][0];
    *pred1 += nn * ys_mean[j][1];
  }
  if (count != 0) {
    *pred0 /= count;
    *pred1 /= count;
  }
}


void mean_tx(double *y, int *tx, double *w, int *leaf, int ell, int *n, double *pred0, double *pred1) {
  double s_count[2] = {0, 0};

  *pred0 = 0;
  *pred1 = 0;

  for (int i = 0; i < *n; ++i) {
    if (leaf[i] == ell) {
      s_count[tx[i]] += w[i];
      if (tx[i] == 0) {
        *pred0 += w[i] * y[i];
      }
      else {
        *pred1 += w[i] * y[i];
      }
    }
  }

  *pred0 /= s_count[0];
  *pred1 /= s_count[1];

}

// no honest
struct Tree {
  int numNodes;
  int *var;
  double *val;
  int *left;
  int *right;
  int *leaf;
  double *pred;
  double *cost;
};

struct Tree createTreeStruct(int numNodes) {
  struct Tree tree;
  tree.numNodes = numNodes;
  tree.var = (int *)malloc(sizeof(int) * numNodes);
  tree.val = (double *)malloc(sizeof(double) * numNodes);
  tree.left = (int *)malloc(sizeof(int) * numNodes);
  tree.right = (int *)malloc(sizeof(int) * numNodes);
  tree.leaf = (int *)malloc(sizeof(int) * numNodes);
  tree.pred = (double *)malloc(sizeof(double) * numNodes * 2);
  tree.cost = (double *)malloc(sizeof(double) * numNodes);

  for (int i = 0; i < numNodes; ++i) {
    tree.var[i] = -1;
  }
  return tree;
}


// no error checking of 0/1 assignment, no randomn, no txbar, no randomSplit, no alpha, no randomp, no *ps, no *Xe, *ye, *txe, *pse
// no Xe yet
// void bestSplitTree(double *X, double *y, int *tx, double *Xe, double *ye, int *txe, double *pse,
//                    int *maxLeaves, int *usePropensity, int *s, int *se, int *isConstVar, int *n, int *p, int *ne,
//                    int *var, double *val, int *left, int *right, double *pred0, double *pred1, double *cost) {
void bestSplitTree(double *X, double *y, int *tx, 
                   double *Xe, double *ye, int *txe,
                   int *maxLeaves, int *usePropensity, int *s, int *se, int *isConstVar, int *n, int *p, int *ne,
                   double *splitSpread, 
                   int *var, double *val, int *left, int *right, double *pred0, double *pred1, double *cost,
                   double *pred0e, double *pred1e) {
  int num_s, *ns = &num_s;
  if (*s != -1) {
    *ns = pmax_int(s, *n) + 1;
  }
  else {
    *ns = -1;  // = to R NULL
  }


  double *w = (double *)malloc((*n) * sizeof(double));
  double *we = (double *)malloc((*ne) * sizeof(double));
  for (int i = 0; i < *n; ++i) {
    w[i] = 1;
  }
  for (int i = 0; i < *ne; ++i) {
    we[i] = 1;
  }
  // if (*ps == NA_REAL) {
  //   for (int i = 0; i < *n; ++i) {
  //     w[i] = 1;
  //   }
  // }
  // else {
  //   for (int i = 0; i < *n; ++i) {
  //     w[i] = 1.0 / ps[i] * tx[i] + 1.0 / (1 - ps[i]) * (1 - tx[i]);
  //   }
  // }

  // double we[*ne];
  // if (*Xe != NA_REAL) {
  //   if (*ps == NA_INTEGER) {
  //     for (int i = 0; i < *ne; ++i) {
  //       we[i] = 1;
  //     }
  //   }
  //   else {
  //     for (int i = 0; i < *ne; ++i) {
  //       we[i] = 1.0 / pse[i] * txe[i] + 1.0 / (1 - pse[i]) * (1 - txe[i]);
  //     }
  //   }
  // }

  int *order = (int *)malloc((*p) * (*n) * sizeof(int));
  // int order[(*p) * (*n)];  // here order starts from 0
  for (int j = 0; j < *p; ++j) {
    order_index_double(X + j * (*n), order + j * (*n), *n);
  }

  int kSplits = 0;
  int ell;

  int *leaf = (int *)malloc((*n) * sizeof(int));
  int *leafe = (int *)malloc((*ne) * sizeof(int));
  int *leaf_var = (int *)malloc((*n) * sizeof(int));

  double *leaf_val = (double *)malloc((*n) * sizeof(double));
  double *leaf_criterion = (double *)malloc((*n) * sizeof(double));

  // int leaf[*n], leafe[*ne], leaf_var[*n];
  // double leaf_val[*n], leaf_criterion[*n];
//  double pred[*n * 2];
//  double cost[*n];
//  int var[*n], left[*n], right[*n];
//  double val[*n];
  int nextAvailable = 1, numNodes = 0;

  for (int i = 0; i < *n; ++i) {
    leaf[i] = 0;
    leaf_criterion[i] = -1;
  }
  for (int i = 0; i < *ne; ++i) {
    leafe[i] = 0;
  }

  while (1) {
    int node, variable;
    double value, objective;
    if (kSplits + 1 < *maxLeaves) {
      if (*usePropensity) {
        findBestSplit_propensity(X, y, tx, w, order, n, p, leaf, leaf_var, leaf_val, leaf_criterion, s, ns, isConstVar, 
          &node, &variable, &value, &objective, splitSpread);
      }
      else {
        findBestSplit(X, y, tx, w, order, n, p, leaf, leaf_var, leaf_val, leaf_criterion, &node, &variable, &value, &objective, splitSpread);
      }
      ell = node;
    }
    else {
      ell = -1;  // NA in Rcode
    }


    if (ell == -1 || ell == NA_INTEGER) {
      int *uniqueLeaves = (int *)malloc((*n) * sizeof(int));
      // int uniqueLeaves[*n];
      for (int i = 0; i < *n; ++i) {
        uniqueLeaves[i] = leaf[i];
      }
      int nu = findUnique(uniqueLeaves, *n);

      for (int k = 0; k < nu; ++k) {
        ell = uniqueLeaves[k];
        if (*usePropensity) {
          mean_propensity_tx(y, tx, w, leaf, ell, n, s, ns, pred0 + ell, pred1 + ell);
        }
        else {
          mean_tx(y, tx, w, leaf, ell, n, pred0 + ell, pred1 + ell);
        }

        if (*txe != -1) {
          if (*usePropensity) {
            mean_propensity_tx(ye, txe, we, leafe, ell, ne, se, ns, pred0e + ell, pred1e + ell);
          }
          else {
            mean_tx(ye, txe, we, leafe, ell, ne, pred0e + ell, pred1e + ell);
          }
        }

        cost[ell] = 0;
        for (int i = 0; i < *n; ++i) {
          if (leaf[i] == ell) {
            if (tx[i] == 0) {
              cost[ell] += w[i] * (y[i] - pred0[ell]) * (y[i] - pred0[ell]);
            }
            else {
              cost[ell] += w[i] * (y[i] - pred1[ell]) * (y[i] - pred1[ell]);
            }
          }
        }
      }
      free(uniqueLeaves);
      break;
    }


    kSplits++;

    if (*usePropensity) {
      mean_propensity_tx(y, tx, w, leaf, ell, n, s, ns, pred0 + ell, pred1 + ell);
    }
    else {
      mean_tx(y, tx, w, leaf, ell, n, pred0 + ell, pred1 + ell);
    }

    if (*txe != -1) {
      if (*usePropensity) {
        mean_propensity_tx(ye, txe, we, leafe, ell, ne, se, ns, pred0e + ell, pred1e + ell);
      }
      else {
        mean_tx(ye, txe, we, leafe, ell, ne, pred0e + ell, pred1e + ell);
      }
    }

    cost[ell] = 0;
    for (int i = 0; i < *n; ++i) {
      if (leaf[i] == ell) {
        if (tx[i] == 0) {
          cost[ell] += w[i] * (y[i] - pred0[ell]) * (y[i] - pred0[ell]);
        }
        else {
          cost[ell] += w[i] * (y[i] - pred1[ell]) * (y[i] - pred1[ell]);
        }
      }
    }

    var[ell] = variable;
    val[ell] = value;

    left[ell] = nextAvailable;
    right[ell] = nextAvailable + 1;
    numNodes = nextAvailable + 2;
    nextAvailable += 2;
    for (int i = 0; i < *n; ++i) {
      if (leaf[i] == ell) {
        if (X[var[ell] * (*n) + i] < val[ell]) {
          leaf[i] = left[ell];
        }
        else {
          leaf[i] = left[ell] + 1;
        }
      }
    }
    if (*txe != -1) {
      for (int i = 0; i < *ne; ++i) {
        if (leafe[i] == ell) {
          if (Xe[var[ell] * (*ne) + i] < val[ell]) {
            leafe[i] = left[ell];
          }
          else {
            leafe[i] = left[ell] + 1;
          }
        }
      }
    }
  }

  // struct Tree tree = createTreeStruct(numNodes);

  // for (int j = 0; j < numNodes; ++j) {
  //   tree.var[j] = var[j];
  //   tree.val[j] = val[j];
  //   tree.left[j] = left[j];
  //   tree.right[j] = right[j];
  //   tree.cost[j] = cost[j];
  //   for (int k = 0; k < 2; ++k) {
  //     tree.pred[j * 2] = pred[j * 2];
  //     tree.pred[j * 2 + 1] = pred[j * 2 + 1];
  //   }
  // }
  // return tree;

  free(w);
  free(we);
  free(order);
  free(leaf);
  free(leafe);
  free(leaf_var);
  free(leaf_val);
  free(leaf_criterion);
}



void predictTree(double *newX, int *n, double *result0, double *result1,
                 int *var, double *val, int *left, int *right, double *pred0, double *pred1) {
  for (int i = 0; i < *n; ++i) {
    int j = 0;
    while (var[j] != NA_INTEGER) {
      if (newX[var[j] * (*n) + i] >= val[j]) {
        j = right[j];
      }
      else {
        j = left[j];
      }
    }
    result0[i] = pred0[j];
    result1[i] = pred1[j];
  }
}




void initializeResults(int nTrees, int maxLeaves, int nRow, int *var, double *val, int *left, int *right, 
                       double *pred0, double *pred1, double *cost, double *pred0e, double *pred1e) {
  int maxNodes = 2 * maxLeaves - 1;
  for (int i = 0; i < nTrees * maxNodes; ++i) {
    var[i] = NA_INTEGER;
    val[i] = NA_REAL;
    left[i] = NA_INTEGER;
    right[i] = NA_INTEGER;
    pred0[i] = NA_REAL;
    pred1[i] = NA_REAL;
    cost[i] = NA_REAL;
    pred0e[i] = NA_REAL;
    pred1e[i] = NA_REAL;
  }
}

void causalBoosting(double *X, double *y, int *tx, 
                     double *Xe, double *ye, int *txe,
                     int *nTrees, int *maxLeaves, double *eps,
                     int *usePropensity, int *s, int *se, int *isConstVar, 
                     int *n, int *p, int *ne, double *vtxeff, 
                     double *splitSpread,
                     int *var, double *val, int *left, int *right, double *y0bar, double *y1bar, double *pred0, double *pred1, double *cost,
                     double *pred0e, double *pred1e, 
                     double *G0, double *G1, double *err_y, double *err, double *tauhat) {

  initializeResults(*nTrees, *maxLeaves, *n, var, val, left, right, pred0, pred1, cost, pred0e, pred1e);

  int maxNodes = 2 * (*maxLeaves) - 1;

  double *r = (double *)malloc((*n) * sizeof(double));
  double *re = (double *)malloc((*ne) * sizeof(double));

  // *y0bar = meanByFactor(y, tx, 0, *n);
  // *y1bar = meanByFactor(y, tx, 1, *n);
  
  // for (int i = 0; i < *n; ++i) {
  //   if (tx[i] == 0) {
  //     r[i] = y[i] - *y0bar;
  //   }
  //   else {
  //     r[i] = y[i] - *y1bar;
  //   }
  //   G0[i] = *y0bar;
  //   G1[i] = *y1bar;
  // }
  
  for (int i = 0; i < *n; ++i) {
    r[i] = y[i];
  }
  for (int i = 0; i < *ne; ++i) {
    re[i] = ye[i];
  }

  // for (int k = 0; k < *nTrees; ++k) {
  //   bestSplitTree(X, r, tx, Xe, ye, txe, maxLeaves, usePropensity, s, se, isConstVar, n, p, ne,
  //                 splitSpread, 
  //                 var + k * maxNodes, val + k * maxNodes, left + k * maxNodes, right + k * maxNodes, 
  //                 pred0 + k * maxNodes, pred1 + k * maxNodes, cost + k * maxNodes,
  //                 pred0e + k * maxNodes, pred1e + k * maxNodes);
  //   double fitted0[*n], fitted1[*n];
  //   predictTree(X, n, fitted0, fitted1, 
  //               var + k * maxNodes, val + k * maxNodes, left + k * maxNodes, right + k * maxNodes, 
  //               pred0 + k * maxNodes, pred1 + k * maxNodes);
  //   err[k] = err_y[k] = 0;
  //   for (int i = 0; i < *n; ++i) {
  //     G0[i] += *eps * fitted0[i];
  //     G1[i] += *eps * fitted1[i];
  //     if (tx[i] == 0) {
  //       r[i] -= *eps * fitted0[i];
  //       err_y[k] += (y[i] - G0[i]) * (y[i] - G0[i]);
  //     }
  //     else {
  //       r[i] -= *eps * fitted1[i];
  //       err_y[k] += (y[i] - G1[i]) * (y[i] - G1[i]);
  //     }
  //     tauhat[*n * k + i] = G1[i] - G0[i];
  //     err[k] += (tauhat[*n * k + i] - *vtxeff) * (tauhat[*n * k + i] - *vtxeff);
  //   }
  //   err[k] /= *n;
  //   err_y[k] /= *n;
  // }

  for (int k = 0; k < *nTrees; ++k) {
    bestSplitTree(X, r, tx, Xe, re, txe, maxLeaves, usePropensity, s, se, isConstVar, n, p, ne,
                  splitSpread, 
                  var + k * maxNodes, val + k * maxNodes, left + k * maxNodes, right + k * maxNodes, 
                  pred0 + k * maxNodes, pred1 + k * maxNodes, cost + k * maxNodes,
                  pred0e + k * maxNodes, pred1e + k * maxNodes);


    double *fitted0 = (double *)malloc((*n) * sizeof(double));
    double *fitted1 = (double *)malloc((*n) * sizeof(double));
    double *fitted0e = (double *)malloc((*ne) * sizeof(double));
    double *fitted1e = (double *)malloc((*ne) * sizeof(double));


    predictTree(X, n, fitted0, fitted1, 
                var + k * maxNodes, val + k * maxNodes, left + k * maxNodes, right + k * maxNodes, 
                pred0 + k * maxNodes, pred1 + k * maxNodes);

    if (*txe != -1) {
      predictTree(Xe, ne, fitted0e, fitted1e,
                  var + k * maxNodes, val + k * maxNodes, left + k * maxNodes, right + k * maxNodes, 
                  pred0e + k * maxNodes, pred1e + k * maxNodes);
    }

    err[k] = err_y[k] = 0;
    for (int i = 0; i < *n; ++i) {
      G0[i] += *eps * fitted0[i];
      G1[i] += *eps * fitted1[i];
      if (tx[i] == 0) {
        r[i] -= *eps * fitted0[i];
        err_y[k] += (y[i] - G0[i]) * (y[i] - G0[i]);
      }
      else {
        r[i] -= *eps * fitted1[i];
        err_y[k] += (y[i] - G1[i]) * (y[i] - G1[i]);
      }
      tauhat[*n * k + i] = G1[i] - G0[i];
      err[k] += (tauhat[*n * k + i] - *vtxeff) * (tauhat[*n * k + i] - *vtxeff);
    }

    err[k] /= *n;
    err_y[k] /= *n;

    if (*txe != -1) {
      for (int i = 0; i < *ne; ++i) {
        // G0[i] += *eps * fitted0[i];
        // G1[i] += *eps * fitted1[i];
        if (txe[i] == 0) {
          re[i] -= *eps * fitted0e[i];
        }
        else {
          re[i] -= *eps * fitted1e[i];
        }
        // tauhat[*n * k + i] = G1[i] - G0[i];
        // err[k] += (tauhat[*n * k + i] - *vtxeff) * (tauhat[*n * k + i] - *vtxeff);
      }
    }
    free(fitted0);
    free(fitted1);
    free(fitted0e);
    free(fitted1e);
  }
  free(r);
  free(re);
}
