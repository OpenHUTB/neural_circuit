double
fmod(double da, double db)
{
	double n;

	n = da/db;
	return(da - n*db);
}
