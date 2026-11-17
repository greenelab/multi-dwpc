import { useEffect, useId, useRef, useState } from "react";

interface ComboboxProps<T> {
  label: string;
  value: string;
  options: T[];
  getKey: (option: T) => string;
  getLabel: (option: T) => string;
  onSelect: (option: T) => void;
  /** Set for server-side search; otherwise options are filtered locally. */
  onSearch?: (text: string) => void;
  loading?: boolean;
  placeholder?: string;
}

/** Searchable select: the current value shows until the user starts typing. */
export function Combobox<T>({
  label,
  value,
  options,
  getKey,
  getLabel,
  onSelect,
  onSearch,
  loading = false,
  placeholder = "Search...",
}: ComboboxProps<T>) {
  const id = useId();
  const inputRef = useRef<HTMLInputElement>(null);
  const listRef = useRef<HTMLUListElement>(null);
  const [open, setOpen] = useState(false);
  const [text, setText] = useState("");
  const [activeIndex, setActiveIndex] = useState(0);

  const needle = text.trim().toLowerCase();
  const visible =
    onSearch || !needle
      ? options
      : options.filter((option) => getLabel(option).toLowerCase().includes(needle));

  useEffect(() => {
    if (open) {
      listRef.current?.children[activeIndex]?.scrollIntoView({ block: "nearest" });
    }
  }, [open, activeIndex]);

  const updateText = (next: string) => {
    setText(next);
    setActiveIndex(0);
    onSearch?.(next);
  };

  const close = () => {
    setOpen(false);
    updateText("");
  };

  const choose = (option: T) => {
    onSelect(option);
    close();
    inputRef.current?.blur();
  };

  const onKeyDown = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === "ArrowDown" || event.key === "ArrowUp") {
      event.preventDefault();
      setOpen(true);
      const step = event.key === "ArrowDown" ? 1 : -1;
      setActiveIndex((i) => Math.min(Math.max(i + step, 0), Math.max(visible.length - 1, 0)));
    } else if (event.key === "Enter" && open && visible[activeIndex] !== undefined) {
      event.preventDefault();
      choose(visible[activeIndex]);
    } else if (event.key === "Escape") {
      close();
      inputRef.current?.blur();
    }
  };

  return (
    <div className="field combobox">
      <label htmlFor={id}>{label}</label>
      <div className="combobox-control">
        <input
          ref={inputRef}
          id={id}
          role="combobox"
          aria-expanded={open}
          aria-controls={`${id}-list`}
          aria-autocomplete="list"
          aria-activedescendant={open && visible.length ? `${id}-${activeIndex}` : undefined}
          autoComplete="off"
          title={value}
          value={open ? text : value}
          placeholder={value || placeholder}
          onChange={(event) => updateText(event.target.value)}
          onFocus={() => setOpen(true)}
          onClick={() => setOpen(true)}
          onBlur={close}
          onKeyDown={onKeyDown}
        />
        <svg className="combobox-chevron" viewBox="0 0 24 24" aria-hidden="true">
          <path d="M7.4 8.6 12 13.2l4.6-4.6L18 10l-6 6-6-6z" />
        </svg>
      </div>
      {open && (
        <ul ref={listRef} id={`${id}-list`} role="listbox" aria-label={label}>
          {visible.map((option, index) => (
            <li
              key={getKey(option)}
              id={`${id}-${index}`}
              role="option"
              aria-selected={index === activeIndex}
              onMouseEnter={() => setActiveIndex(index)}
              onMouseDown={(event) => {
                event.preventDefault(); // keep input focus until the choice lands
                choose(option);
              }}
            >
              {getLabel(option)}
            </li>
          ))}
          {!visible.length && <li className="combobox-empty">{loading ? "Searching..." : "No results"}</li>}
        </ul>
      )}
    </div>
  );
}
