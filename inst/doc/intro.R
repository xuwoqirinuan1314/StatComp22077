## ----eval=FALSE---------------------------------------------------------------
#  hessian_i = function(
#    data_i, #i^th observation: (Y_i, X_i, Z_i)
#    symbolic_hessian, #matrix with elements that are the symbolic expression for the derivatives
#    param_vals, #estimates of first stage and second stage parameter estimates (theta and beta)
#    dimension #the dimension of the hessian matrix
#    ){
#  
#    env = as.list(c(data_i, param_vals))
#  
#    # Iterate through every element in symbolic_hessian and evaluate it
#    # at the i^th observation using the given parameter estimates
#    hessian = matrix(
#      sapply(symbolic_hessian, eval, env = env, enclos = NULL),
#      dimension)
#  
#    return(hessian)
#  
#  }

## ----eval=FALSE---------------------------------------------------------------
#  int part(int* r, int low, int hight)
#  {
#    int i = low, j = hight, pivot = r[low];
#    while (i < j)
#    {
#      while (i<j && r[j]>pivot)
#      {
#        j--;
#      }
#      if (i < j)
#      {
#        swap(r[i++], r[j]);
#      }
#      while (i < j && r[i] <= pivot)
#      {
#        i++;
#      }
#      if (i < j)
#      {
#        swap(r[i], r[j--]);
#      }
#    }
#    return i;
#  }
#  void Quicksort(int* r, int low, int hight)
#  {
#    int mid;
#    if (low < hight)
#    {
#      mid = part(r, low, hight);
#      Quicksort(r, low, mid - 1);
#      Quicksort(r, mid+1, hight);
#    }
#  }
#  int main()
#  {
#    int a[10001];
#    int  N;
#    cout << "Please enter the number of data to sort: " << endl;
#    cin >> N;
#    cout << "Please enter the data to sort: " << endl;
#    for (int i = 0; i < N; i++)
#    {
#      cin >> a[i];
#    }
#    cout << endl;
#    Quicksort(a, 0, N - 1);
#    cout << "The sorted sequence is: " << endl;
#    for (int i = 0; i < N; i++)
#    {
#      cout << a[i] << " ";
#    }
#    cout << endl;
#    return 0;
#  }

