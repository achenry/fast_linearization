linturb_desc = FASTLinearizationFile()

x = {k: v for (k, v) in zip(x_desc, x_op) if k[0:2] == 'ED'}
x_desc = linturb_desc.short_descr(linturb_models.DescStates)
u = {k: v for (k, v) in zip(u_desc, u_op) if k[0:2] == 'ED'}
u_desc = linturb_desc.short_descr(linturb_models.DescCntrlInpt)

y = {k: v for (k, v) in zip(y_desc, y_op) if k[0:2] == 'ED'}
y_desc = linturb_desc.short_descr(linturb_models.DescOutput)


