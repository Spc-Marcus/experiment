
def add_WarmStart(model, rows, cols):
        for row in rows:
            if row in model.lp_rows:
                model.lp_rows[row][0].start = 1
        for col in cols:
            if col in model.lp_cols:
                model.lp_cols[col][0].start = 1
        model.update()
        return model